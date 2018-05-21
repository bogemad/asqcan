#!/usr/bin/env python2

import sys, os, subprocess, multiprocessing, logging, shutil, gzip, re
from StringIO import StringIO
from collections import defaultdict

class Asqcan():
	def __init__(self, threads, mem, db, reads_dir, outdir, verbosity):
		self.threads = str(threads)
		self.mem = str(mem)
		self.db = db
		self.reads_dir = reads_dir
		self.reads_list = [ os.path.join(reads_dir, reads) for reads in os.listdir(reads_dir) if (os.path.isfile(os.path.join(reads_dir, reads)) and (reads.endswith('.fastq') or reads.endswith('.fastq.gz')))]
		self.outdir = outdir
		self.fastqc_dir = os.path.join(outdir, 'qc_plots', 'fastqc')
		self.raw_ass_dir = os.path.join(outdir, 'raw_assemblies')
		self.quast_dir = os.path.join(outdir, 'qc_plots', 'quast')
		self.blobtools_dir = os.path.join(outdir, 'qc_plots', 'blobtools')
		self.ann_ass_dir = os.path.join(outdir, 'annotated_assemblies')
		self.temp_dir = os.path.join(outdir, '.temp')
		self.successful = defaultdict(list)
		self.failed = defaultdict(list)
		self.verbosity = verbosity
		self.steps = [
					'fastqc',
					'spades',
					'quast',
					'blastn',
					'blobtools create',
					'blobtools view',
					'blobtools plot genus',
					'blobtools plot species',
					'prokka'
					]

	
	def open_reads(self, reads):
		if reads.endswith('.fastq'):
			reads_handle = open(reads, 'r')
		elif reads.endswith('fastq.gz'):
			reads_handle = gzip.open(reads, 'rt')
		else:
			logging.info("Reads: {} not fastq!")
			sys.exit(1)
		return reads_handle
	
	def clear_temp_dir(self):
		shutil.rmtree(self.temp_dir)
		os.makedirs(self.temp_dir)
	
	def kmer_calc(self, reads):
		sample_name = get_sample_name(reads)
		reads_handle = self.open_reads(reads)
		max_read_len = max(len(line) for i, line in enumerate(reads_handle) if (i % 4) == 1 )
		reads_handle.close()
		run_ks = []
		for k in [21,33,55,77,99,127]:
			run_ks.append(str(k))
			if k > (max_read_len/2):
				logging.info("{}: Running spades with kmers {}".format(sample_name, ",".join(run_ks)))
				return ','.join(run_ks)
	
	def check_fix_SRA_reads(self, reads):
		sample_name = get_sample_name(reads)
		reads_handle = self.open_reads(reads)
		headers = []
		for i, l in enumerate(reads_handle):
			if i % 4 == 0:
				headers.append(l.strip().split()[0])
			if i == 7:
				break
		reads_handle.close()
		if headers[0] != headers[1]:
			sO = re.match(r'^@[SED]RR[0-9]+', headers[0])
			if sO:
				logging.info("{}: Incompatible read archive format detected. Attempting to fix...".format(sample_name))
				make_outdirs(self.temp_dir)
				reads_handle = self.open_reads(reads)
				fixed_reads = os.path.join(self.temp_dir, "{}.fixed.fastq".format(sample_name))
				with open(fixed_reads, 'w') as fixed_reads_handle:
					for i, line in enumerate(reads_handle):
						if i % 4 == 0 or i % 4 == 2:
							line = re.sub(r'^((@|\+)[SED]RR[^.]+\.[^.]+)\.(1|2)', r'\1 \3', line)
						fixed_reads_handle.write(line)
				logging.info("{}: Reads fix complete.".format(sample_name))
				return fixed_reads
			else:
				logging.info("{}: Incompatible read header format detected. Headers (writing between @ and first space) for paired reads do not match. This will caused problems for bwa in the spades pipeline when run with the --careful option. If you wish to use this option for assemblies, please correct your fastq file headers.".format(sample_name))
		return reads
	
	def assemble_reads(self, reads):
		sample_name = get_sample_name(reads)
		make_outdirs(self.raw_ass_dir)
		output = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
		os.path.isfile(output)
		if not os.path.isfile(output):
			reads = self.check_fix_SRA_reads(reads)
			logging.info("{}: Assembling reads with spades...".format(sample_name))
			if os.path.isdir(os.path.join(self.raw_ass_dir, sample_name)):
				shutil.rmtree(os.path.join(self.raw_ass_dir, sample_name))
			exitcode = run_command([
									'spades.py',
									'--careful',
									'--pe1-12', reads,
									'-t', self.threads,
									'-m', self.mem,
									'-k', self.kmer_calc(reads),
									'--cov-cutoff', 'auto',
									'-o', os.path.join(self.raw_ass_dir, sample_name)
									])
			if exitcode == 1:
				logging.info("{}: Retrying spades assembly without --careful option".format(sample_name))
				if os.path.isdir(os.path.join(self.raw_ass_dir, sample_name)):
					shutil.rmtree(os.path.join(self.raw_ass_dir, sample_name))
				exitcode = run_command([
										'spades.py',
										'--pe1-12', reads,
										'-t', self.threads,
										'-m', self.mem,
										'-k', self.kmer_calc(reads),
										'--cov-cutoff', 'auto',
										'-o', os.path.join(self.raw_ass_dir, sample_name)
										])
			self.log_exit_status(exitcode, 'spades', sample_name)
			self.clear_temp_dir()
	
	def run_quast(self, reads):
		sample_name = get_sample_name(reads)
		make_outdirs(self.quast_dir)
		output = os.path.join(self.quast_dir, sample_name, 'report.pdf')
		if not os.path.isfile(output):
			logging.info("{}: Quality checking assembly with quast...".format(sample_name))
			available_dirs = [ dir for dir in ['K21','K33','K55','K77','K99','K127'] if os.path.isdir(os.path.join(self.raw_ass_dir, sample_name, dir))]
			ass_list = [ os.path.join(self.raw_ass_dir, sample_name, dir, 'final_contigs.fasta') for dir in available_dirs if 'final_contigs.fasta' in os.listdir(os.path.join(self.raw_ass_dir, sample_name, dir)) ] + [os.path.join(self.raw_ass_dir, sample_name, 'contigs.fasta'), os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')]
			exitcode = run_command(([
									'quast.py',
									'-t', self.threads,
									'-o', os.path.join(self.quast_dir, sample_name)
									] + ass_list))
			self.log_exit_status(exitcode, 'quast', sample_name)
	
	def run_fastqc(self, reads):
		make_outdirs(self.fastqc_dir)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.fastqc_dir, "{}_fastqc.html".format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without fastqc reports. Starting fastqc...".format(len(not_done))) 
			exitcode = run_command([
									'parallel', '-j{}'.format(self.threads),
									'fastqc',
									'--outdir={}'.format(self.fastqc_dir),
									'{}', ':::'] + not_done)
			if exitcode == 0:
				for reads in not_done:
					sample_name = get_sample_name(reads)
					logging.info("{}: Execution of fastqc completed successfully.".format(sample_name))
					self.successful['fastqc'].append(sample_name)
			else:
				for reads in not_done:
					sample_name = get_sample_name(reads)
					logging.info("{}: Execution of fastqc failed.".format(sample_name))
					if self.verbosity == False:
						logging.info("For more information please rerun the pipeline with the -f|--force and -v|--verbose options.")
					self.failed['fastqc'].append(sample_name)
			logging.info("\nExecution of fastqc completed. Result summary:\n\nSuccessful: {}\nFailed: {}\n\nBeginning assembly of reads...".format(len(self.successful['fastqc']), len(self.failed['fastqc'])))
		elif len(not_done) == 0 and len(self.reads_list) > 0:
			for reads in self.reads_list:
				sample_name = get_sample_name(reads)
				self.successful['fastqc'].append(sample_name)
		else:
			logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
			sys.exit(1)
	
	def run_blastn(self, reads):
		sample_name = get_sample_name(reads)
		make_outdirs(os.path.join(self.blobtools_dir, sample_name))
		assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
		blast_hits_file = os.path.join(self.blobtools_dir, sample_name, '{}.blast_hits'.format(sample_name))
		if not os.path.isfile(blast_hits_file):
			logging.info("{}: Quality checking assembly with blobtools...".format(sample_name))
			logging.info("{}: blobtools: Searching nt database with blastn...".format(sample_name))
			exitcode = run_command([
									'blastn', 
									'-query', assembly,
									'-db', self.db,
									'-outfmt', '6 qseqid staxids bitscore std',
									'-max_target_seqs', '10',
									'-max_hsps', '1',
									'-evalue', '1e-25',
									'-num_threads', self.threads,
									'-out', blast_hits_file
									])
			self.log_exit_status(exitcode, 'blastn', sample_name)
	
	def run_bt_create(self, reads):
		sample_name = get_sample_name(reads)
		assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
		blast_hits_file = os.path.join(self.blobtools_dir, sample_name, '{}.blast_hits'.format(sample_name))
		blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
		if not os.path.isfile(blobdb_file):
			logging.info("{}: blobtools: Building blobDB...".format(sample_name))
			logging.info("Generating {} blobDB...".format(sample_name))
			exitcode = run_command([
						'blobtools', 'create',
									'-i', assembly, 
									'-y', 'spades',
									'-t', blast_hits_file,
									'-x', 'bestsumorder',
									'-o', os.path.join(self.blobtools_dir, sample_name, sample_name)
									])
			self.log_exit_status(exitcode, 'blobtools create', sample_name)
	
	def run_bt_view(self, reads):
		sample_name = get_sample_name(reads)
		blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
		blob_table = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.bestsumorder.table.txt'.format(sample_name))
		if not os.path.isfile(blob_table):
			logging.info("{}: blobtools: Generating blob table...".format(sample_name))
			exitcode = run_command([
									'blobtools', 'view',
									'-i', blobdb_file,
									'-x', 'bestsumorder',
									'-r', 'species',
									'-b',
									'-o', os.path.join(self.blobtools_dir, sample_name, sample_name)
									])
			self.log_exit_status(exitcode, 'blobtools view', sample_name)
	
	def run_bt_plot_genus(self, reads):
		sample_name = get_sample_name(reads)
		blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
		blob_plot_genus = os.path.join(self.blobtools_dir, sample_name, '{}.bestsumorder.genus.p7.span.100.blobplot.spades.png'.format(blobdb_file))
		if not os.path.isfile(blob_plot_genus):
			logging.info("{}: blobtools: Generating blob genus plot...".format(sample_name))
			exitcode = run_command([
									'blobtools', 'plot',
									'-i', blobdb_file,
									'-x', 'bestsumorder',
									'-r', 'genus',
									'-l', '100',
									'-o', os.path.join(self.blobtools_dir, sample_name, sample_name)
									])
			self.log_exit_status(exitcode, 'blobtools plot genus', sample_name)
	
	def run_bt_plot_species(self, reads):
		sample_name = get_sample_name(reads)
		blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
		blob_plot_species = os.path.join(self.blobtools_dir, sample_name, '{}.bestsumorder.species.p7.span.100.blobplot.spades.png'.format(blobdb_file))
		if not os.path.isfile(blob_plot_species):
			logging.info("{}: blobtools: Generating blob species plot...".format(sample_name))
			exitcode = run_command([
									'blobtools', 'plot',
									'-i', blobdb_file,
									'-x', 'bestsumorder',
									'-r', 'species',
									'-l', '100',
									'-o', os.path.join(self.blobtools_dir, sample_name, sample_name)
									])
			self.log_exit_status(exitcode, 'blobtools plot species', sample_name)
	
	def run_prokka(self, reads):
		sample_name = get_sample_name(reads)
		make_outdirs(self.ann_ass_dir)
		assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
		output =  os.path.join(self.ann_ass_dir, sample_name, '{}.gbk'.format(sample_name))
		if not os.path.isfile(output):
			logging.info("{}: Annotating assembly with prokka...".format(sample_name))
			exitcode = run_command([
									'prokka',
									'--outdir', os.path.join(self.ann_ass_dir, sample_name),
									'--prefix', sample_name,
									'--cpus', self.threads,
									assembly
									])
			self.log_exit_status(exitcode, 'prokka', sample_name)
	
	def generate_report(self, report_path):
		out_d = defaultdict(list)
		for reads in self.reads_list:
			sample_name = get_sample_name(reads)
			out_d[sample_name] = []
			for step in self.steps:
				if sample_name in self.successful[step]:
					out_d[sample_name].append("Completed")
				elif sample_name in self.failed[step]:
					out_d[sample_name].append("Failed")
				else:
					out_d[sample_name].append("Not_run")
		with open(report_path, 'w') as report_handle:
			report_handle.write('\t'.join(['sample_name'] + self.steps) + '\n')
			for key, val in sorted(out_d.items()):
				report_handle.write('\t'.join([key] + val) + '\n')
	
	def import_previous_report(self, report_path):
		logging.info("Previous report found. Importing...")
		with open(report_path) as report_handle:
			for i, line in enumerate(report_handle):
				if i == 0:
					continue
				line_l = line.strip().split()
				header_l = ['sample_name'] + self.steps
				for j, status in enumerate(line_l):
					if j == 0:
						continue
					if status == "Completed":
						self.successful[header_l[j]].append(line_l[0])
					elif status == "Failed":
						self.failed[header_l[j]].append(line_l[0])
					elif status == "Not_run":
						continue
					else:
						logging.info("Error: Invalid option when importing previous report. Please ensure all items are either 'Completed, Failed, or Not_run")
	
	def run_asqcan(self):
		report_path = os.path.join(self.outdir, "asqcan_report.tsv")
		if os.path.isfile(report_path):
			self.import_previous_report(report_path)
		self.run_fastqc(self.reads_list)
		for reads in self.reads_list:
			self.assemble_reads(reads)
		for reads in self.reads_list:
			if get_sample_name(reads) in self.successful['spades']:
				self.run_quast(reads)
		for reads in self.reads_list:
			if get_sample_name(reads) in self.successful['spades']:
				self.run_blastn(reads)
				self.run_bt_create(reads)
				self.run_bt_view(reads)
				self.run_bt_plot_genus(reads)
				self.run_bt_plot_species(reads)
		for reads in self.reads_list:
			if get_sample_name(reads) in self.successful['spades']:
				self.run_prokka(reads)
		logging.info("Generating summary report, find this at {}".format(os.path.join(self.outdir, "asqcan_report.tsv")))
		self.generate_report(report_path)

	
	def log_exit_status(self, exitcode, step, sample_name):
		if exitcode == 0:
			logging.info("{}: Execution of {} completed successfully.\n".format(sample_name, step))
			self.successful[step].append(sample_name)
		else:
			logging.info("{}: Execution of {} failed.\n".format(sample_name, step))
			if self.verbosity == False:
				logging.info("For more information please rerun the pipeline with the -f|--force and -v|--verbose options.")
			self.failed[step].append(sample_name)
	


def get_sample_name(reads):
	name, suffix = os.path.splitext(os.path.basename(reads))
	if suffix == '.gz':
		name, suffix = os.path.splitext(name)
	return name


def log_subprocess_output(pipe):
	for line in iter(pipe.readline, b''): # b'\n'-separated lines
		logging.debug(line.strip())


def make_outdirs(dir):
	if os.path.isdir(dir) == False:
		os.makedirs(dir)


def run_command(cmd):
	logging.info("Executing command: {}".format(" ".join(cmd)))
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	with process.stdout:
		log_subprocess_output(process.stdout)
	exitcode = process.wait()
	return exitcode


def configure_logs(verbosity, outdir):
	logging.basicConfig(filename=os.path.join(outdir, 'asqcan.log'),format='%(levelname)s:%(message)s',level=logging.DEBUG)
	cmd_log = logging.StreamHandler()
	if verbosity == True:
		cmd_log.setLevel(logging.DEBUG)
	else:
		cmd_log.setLevel(logging.INFO)
	logging.getLogger().addHandler(cmd_log)


def configure_outdir(outdir, force):
	if os.path.isdir(outdir):
		if force == True:
			shutil.rmtree(outdir)
			os.makedirs(outdir)
	else:
		os.makedirs(outdir)

def run(threads, mem, db, reads_dir, force, outdir, verbosity):
	configure_outdir(outdir, force)
	configure_logs(verbosity, outdir)
	logging.info("Beginning asqcan run...")
	job = Asqcan(threads, mem, db, reads_dir, outdir, verbosity)
	if force == True:
		logging.info("Force option given. Deleting previous asqcan run...")
	logging.info("Log for the current run can be found at {}.".format(os.path.join(outdir, 'logs')))
	job.run_asqcan()
	logging.info("Run complete. Thanks your for using asqcan.")
	job.clear_temp_dir()

