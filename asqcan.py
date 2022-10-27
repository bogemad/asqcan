#!/usr/bin/env python

import sys, os, subprocess, multiprocessing, logging, shutil, gzip, re
from io import StringIO
from collections import defaultdict
from distutils.spawn import find_executable

class Asqcan():
    def __init__(self, threads, mem, db, reads_dir, outdir, verbosity, ion_torrent, assembler, searcher):
        self.threads = str(threads)
        self.spades_threads = int(threads/4) if threads >= 4 else 1
        self.mem = str(mem)
        self.db = db
        self.reads_dir = reads_dir
        self.reads_list = [ os.path.join(reads_dir, reads) for reads in os.listdir(reads_dir) if (os.path.isfile(os.path.join(reads_dir, reads)) and reads.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')))]
        self.outdir = outdir
        self.fastp_dir = os.path.join(outdir, 'qc_plots', 'fastp')
        self.raw_ass_dir = os.path.join(outdir, 'raw_assemblies')
        self.quast_dir = os.path.join(outdir, 'qc_plots', 'quast')
        self.blobtools_dir = os.path.join(outdir, 'qc_plots', 'blobtools')
        self.bam_dir = os.path.join(outdir, 'qc_plots', 'bam')
        self.ann_ass_dir = os.path.join(outdir, 'annotated_assemblies')
        self.final_fasta = os.path.join(outdir, 'collated_assemblies', 'genome_fasta')
        self.final_gbk = os.path.join(outdir, 'collated_assemblies', 'genome_gbk')
        self.protein_fasta = os.path.join(outdir, 'collated_assemblies', 'protein_fasta')
        self.temp_dir = os.path.join(outdir, '.temp')
        self.split_reads_dir = os.path.join(self.temp_dir, 'split_reads')
        self.fastp_reads_dir = os.path.join(outdir, 'fastp_reads')
        self.successful = defaultdict(list)
        self.failed = defaultdict(list)
        self.verbosity = verbosity
        self.ion_torrent = ion_torrent
        self.assembler = assembler
        self.searcher = searcher
        self.steps = [
                    'fastp',
                    'assembly',
                    'blast',
                    'blobtools create',
                    'blobtools view',
                    'blobtools plot genus',
                    'blobtools plot species',
                    'prokka'
                    ]

    
    def open_reads(self, reads):
        if reads.endswith('.gz'):
            reads_handle = gzip.open(reads, 'rt')
        else:
            reads_handle = open(reads, 'r')
        return reads_handle
    
    def clear_temp_dir(self):
        if os.path.isdir(self.temp_dir):
            shutil.rmtree(self.temp_dir)
        os.makedirs(self.temp_dir)
    
    def check_fix_SRA_reads(self, reads):
        sample_name = get_sample_name(reads)
        reads_handle = self.open_reads(reads)
        headers = []
        seqs = []
        for i, l in enumerate(reads_handle):
            if i % 4 == 0:
                headers.append(l.strip().split()[0])
            if i % 4 == 1:
                seqs.append(l.strip())
            if i == 7:
                break
        reads_handle.close()
        for seq in seqs:
            if len(seq) > 1000:
                logging.info("{}: Long reads detected, please ensure your reads are sourced from the Illumina sequencing platform.".format(sample_name))
                return False
        if self.ion_torrent == False:
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
                    logging.info("{}: Incompatible read header format detected. Headers (writing between @ and first space) for paired reads do not match. This will caused problems for bwa in the spades pipeline when run with the --careful option. If you wish to use this option for assemblies, please ensure you are using Illumina paired-end reads and correct your fastq file headers.".format(sample_name))
                    return False
        return reads
    
    def assemble_reads_spades(self):
        make_outdirs(self.raw_ass_dir)
        cmds = []
        threads = self.threads
        for reads in self.reads_list:
            sample_name = get_sample_name(reads)
            output = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
            if not os.path.isfile(output):
                reads = self.check_fix_SRA_reads(reads)
                if reads == False:
                    logging.info("{}: Assembly aborted, reads are unsuitable for assembly.".format(sample_name))
                    self.log_exit_status('1', 'assembly', sample_name)
                    return
                if self.ion_torrent == True:
                    clean_reads = os.path.join(self.fastp_reads_dir, '{}.fastq.gz'.format(sample_name))
                    it_opt = ' --iontorrent -s {} '.format(clean_reads)
                else:
                    clean_reads_1 = os.path.join(self.fastp_reads_dir, '{}_1.fastq.gz'.format(sample_name))
                    clean_reads_2 = os.path.join(self.fastp_reads_dir, '{}_2.fastq.gz'.format(sample_name))
                    clean_reads_u = os.path.join(self.fastp_reads_dir, '{}_u.fastq.gz'.format(sample_name))
                    it_opt = ' -1 {} -2 {} -s {} '.format(clean_reads_1, clean_reads_2, clean_reads_u)
                if os.path.isdir(os.path.join(self.raw_ass_dir, sample_name)):
                    shutil.rmtree(os.path.join(self.raw_ass_dir, sample_name))
                cmds.append('spades.py --isolate' + it_opt + ' -t 4 -m {} --cov-cutoff auto -o {}'.format(self.mem, os.path.join(self.raw_ass_dir, sample_name)))
            else:
                exitcode = 0
                self.log_exit_status(exitcode, 'assembly', sample_name)
        if len(cmds) == 0:
            return
        exitcode = run_parallel_commands(cmds, 'spades', str(max(int(self.threads)/4,1)))
        logging.debug(cmds)
        self.clear_temp_dir()
        for reads in self.reads_list:
            if os.path.isfile(os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')):
                self.log_exit_status(exitcode, 'assembly', get_sample_name(reads))
    
    def assemble_reads_skesa(self, reads):
        sample_name = get_sample_name(reads)
        make_outdirs(self.raw_ass_dir)
        output = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
        threads = self.threads
        if not os.path.isfile(output):
            reads = self.check_fix_SRA_reads(reads)
            if reads == False:
                logging.info("{}: Assembly aborted, reads are unsuitable for assembly.".format(sample_name))
                self.log_exit_status('1', 'assembly', sample_name)
                return
            logging.info("{}: Assembling reads with skesa...".format(sample_name))
            if os.path.isdir(os.path.join(self.raw_ass_dir, sample_name)):
                shutil.rmtree(os.path.join(self.raw_ass_dir, sample_name))
            os.mkdir(os.path.join(self.raw_ass_dir, sample_name))
            exitcode = run_command(['skesa', '--cores', threads, '--memory', self.mem, '--fastq', reads, '--use_paired_ends', '--contigs_out', output])
            self.clear_temp_dir()
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'assembly', sample_name)
    
    def run_quast(self):
        make_outdirs(self.quast_dir)
        output = os.path.join(self.quast_dir, 'report.pdf')
        if not os.path.isfile(output):
            logging.info("Quality checking assemblies with quast...")
            assemblies = [ os.path.join(self.raw_ass_dir, dir, 'scaffolds.fasta') for dir in os.listdir(self.raw_ass_dir) ]
            exitcode = run_command(([
                                    'quast',
                                    '-t', self.threads,
                                    '-o', os.path.join(self.quast_dir)
                                    ] + assemblies))
        else:
            exitcode = 0
    
    # def run_fastp(self, reads):
        # make_outdirs(self.fastqc_dir)
        # not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.fastqc_dir, "{}_fastqc.html".format(get_sample_name(reads)))) ]
        # if len(not_done) > 0 and len(self.reads_list) > 0:
            # logging.info("Found {} reads without fastqc reports. Starting fastqc...".format(len(not_done))) 
            # exitcode = run_command([
                                    # 'parallel', '-j{}'.format(self.threads),
                                    # 'fastqc',
                                    # '--outdir={}'.format(self.fastqc_dir),
                                    # '{}', ':::'] + not_done)
            # if exitcode == 0:
                # for reads in not_done:
                    # sample_name = get_sample_name(reads)
                    # logging.info("{}: Execution of fastqc completed successfully.".format(sample_name))
                    # self.successful['fastqc'].append(sample_name)
            # else:
                # for reads in not_done:
                    # sample_name = get_sample_name(reads)
                    # logging.info("{}: Execution of fastqc failed.".format(sample_name))
                    # if self.verbosity == False:
                        # logging.info("For more information please rerun the pipeline with the -f|--force and -v|--verbose options.")
                    # self.failed['fastqc'].append(sample_name)
            # logging.info("\nExecution of fastqc completed. Result summary:\n\nSuccessful: {}\nFailed: {}\n\nBeginning assembly of reads...".format(len(self.successful['fastqc']), len(self.failed['fastqc'])))
        # elif len(not_done) == 0 and len(self.reads_list) > 0:
            # for reads in self.reads_list:
                # sample_name = get_sample_name(reads)
                # self.successful['fastqc'].append(sample_name)
        # else:
            # logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
            # sys.exit(1)
    
    def run_fastp(self):
        make_outdirs(self.fastp_dir)
        make_outdirs(self.fastp_reads_dir)
        not_done = [ get_sample_name(reads) for reads in self.reads_list if not os.path.isfile(os.path.join(self.fastp_dir, "{}_fastp.html".format(get_sample_name(reads)))) ]
        if len(not_done) > 0 and len(self.reads_list) > 0:
            logging.info("Found {} reads without fastp reports. Starting fastp...".format(len(not_done)))
            cmds = []
            if self.ion_torrent == True:
                for inread_file in self.reads_list:
                    sample_name = get_sample_name(inread_file)
                    if not sample_name in not_done:
                        continue
                    outread = os.path.join(self.fastp_reads_dir, '{}.fastq.gz'.format(sample_name))
                    json = os.path.join(self.fastp_dir, "{}_fastp.json".format(sample_name))
                    html = os.path.join(self.fastp_dir, "{}_fastp.html".format(sample_name))
                    cmds.append("fastp -i {0} -o {1} -j {2} -h {3}".format(inread_file, outread, json, html))
            else:
                for inread_file in self.reads_list:
                    sample_name = get_sample_name(inread_file)
                    if not sample_name in not_done:
                        continue
                    json = os.path.join(self.fastp_dir, "{}_fastp.json".format(sample_name))
                    html = os.path.join(self.fastp_dir, "{}_fastp.html".format(sample_name))
                    outread_1 = os.path.join(self.fastp_reads_dir, '{}_1.fastq.gz'.format(sample_name))
                    outread_2 = os.path.join(self.fastp_reads_dir, '{}_2.fastq.gz'.format(sample_name))
                    outread_u = os.path.join(self.fastp_reads_dir, '{}_u.fastq.gz'.format(sample_name))
                    cmds.append("fastp --in1 {0} --interleaved_in --out1 {1} --out2 {2} --unpaired1 {3} --unpaired2 {3} -j {4} -h {5}".format(inread_file, outread_1, outread_2, outread_u, json, html))
            exitcode = run_parallel_commands(cmds, 'fastp', str(max(int(self.threads)/4,1)))
            if exitcode == 0:
                for sample_name in not_done:
                    logging.info("{}: Execution of fastp completed successfully.".format(sample_name))
                    self.successful['fastp'].append(sample_name)
            else:
                for sample_name in not_done:
                    logging.info("{}: Execution of fastp failed.".format(sample_name))
                    if self.verbosity == False:
                        logging.info("For more information please rerun the pipeline with the -f|--force and -v|--verbose options.")
                    self.failed['fastp'].append(sample_name)
            logging.info("\nExecution of fastp completed. Result summary:\n\nSuccessful: {}\nFailed: {}\n\nBeginning assembly of reads...".format(len(self.successful['fastp']), len(self.failed['fastp'])))
        elif len(not_done) == 0 and len(self.reads_list) > 0:
            for reads in self.reads_list:
                sample_name = get_sample_name(reads)
                self.successful['fastp'].append(sample_name)
        else:
            logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
            sys.exit(1)
    
    def run_blastn(self, reads):
        sample_name = get_sample_name(reads)
        make_outdirs(os.path.join(self.blobtools_dir, sample_name))
        assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
        blast_hits_file = os.path.join(self.blobtools_dir, sample_name, '{}.blast_hits'.format(sample_name))
        if not (os.path.isfile(blast_hits_file) or self.db == "Not used"):
            logging.info("{}: Quality checking assembly with blobtools...".format(sample_name))
            logging.info("{}: blobtools: Searching database with blastn...".format(sample_name))
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
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'blast', sample_name)
    
    def run_diamond(self, reads):
        sample_name = get_sample_name(reads)
        make_outdirs(os.path.join(self.blobtools_dir, sample_name))
        assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
        blast_hits_file = os.path.join(self.blobtools_dir, sample_name, '{}.blast_hits'.format(sample_name))
        if not (os.path.isfile(blast_hits_file) or self.db == "Not used"):
            logging.info("{}: Quality checking assembly with blobtools...".format(sample_name))
            logging.info("{}: blobtools: Searching database with diamond...".format(sample_name))
            exitcode = run_command([
                                    'diamond', 'blastx',
                                    '-c1',
                                    '-q', assembly,
                                    '-d', self.db,
                                    '--outfmt', '6', 'qseqid', 'staxids', 'bitscore', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
                                    '--sensitive',
                                    '--max-target-seqs', '1',
                                    '--evalue', '1e-25',
                                    '--threads', self.threads,
                                    '--out', blast_hits_file,
                                    '--taxonmap', os.path.join(os.getenv("HOME"), '.asqcan', 'prot.accession2taxid.gz'),
                                    '--taxonnodes', os.path.join(os.getenv("HOME"), '.asqcan', 'nodes.dmp')
                                    ])
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'blast', sample_name)

        
    def run_bwa_samtools(self, reads):
        sample_name = get_sample_name(reads)
        logging.info("{}: Read mapping to determine assembly coverage...".format(sample_name))
        output = os.path.join(self.bam_dir, '{}.bam'.format(sample_name))
        if not os.path.isfile(output):
            make_outdirs(self.temp_dir)
            logging.info("{}: Building bwa index...".format(sample_name))
            run_command(['bwa', 'index', '-p', os.path.join(self.temp_dir, sample_name), os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')])
            bwa_cmd = ['bwa', 'mem', '-t', self.threads, '-p', os.path.join(self.temp_dir, sample_name), reads ]
            view_cmd = ['samtools', 'view', '-@', str(int(self.threads)-1), '-b' ]
            sort_cmd = ['samtools', 'sort', '-@', str(int(self.threads)-1), '-T', os.path.join(self.temp_dir), '-o', output ]
            logging.info("{}: Generating coverage bam file...".format(sample_name))
            logging.info("Executing command: {}".format("{} | {} | {}".format(' '.join(bwa_cmd), ' '.join(view_cmd), ' '.join(sort_cmd))))
            bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            view_process = subprocess.Popen(view_cmd, stdin=bwa_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sort_process = subprocess.Popen(sort_cmd, stdin=view_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with bwa_process.stderr as bwa_std_err:
                log_subprocess_output(bwa_std_err)
            with view_process.stderr as view_std_err:
                log_subprocess_output(view_std_err)
            with sort_process.stderr as sort_std_err:
                log_subprocess_output(sort_std_err)
            exitcode = sort_process.wait()
            shutil.rmtree(self.temp_dir)
        if not os.path.isfile(output + '.bai'):
            index_cmd = ['samtools', 'index', '-@', str(int(self.threads)-1), output]
            logging.info("Executing command: {}".format("{}".format(' '.join(index_cmd))))
            index_process = subprocess.Popen(index_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            with index_process.stdout as index_std_out:
                log_subprocess_output(index_std_out)
            exitcode2 = index_process.wait()

    
    def run_bt_create(self, reads):
        sample_name = get_sample_name(reads)
        assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
        blast_hits_file = os.path.join(self.blobtools_dir, sample_name, '{}.blast_hits'.format(sample_name))
        blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
        nodes = os.path.join(self.blobtools_dir, 'nodes.dmp')
        names = os.path.join(self.blobtools_dir, 'names.dmp')
        if not (os.path.isfile(nodes) or os.path.isfile(names)):
            run_command(['wget', '-O', os.path.join(self.blobtools_dir, 'taxdmp.zip'), 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip'])
            run_command(['unzip', '-o', '-d', self.blobtools_dir, os.path.join(self.blobtools_dir, 'taxdmp.zip')])
        if not os.path.isfile(blobdb_file):
            logging.info("{}: blobtools: Building blobDB...".format(sample_name))
            logging.info("Generating {} blobDB...".format(sample_name))
            if self.db == "Not used":
                exitcode = run_command([
                                        "blobtools", 'create',
                                        '-i', assembly, 
                                        '-b', os.path.join(self.bam_dir, '{}.bam'.format(sample_name)),
                                        '-x', 'bestsumorder',
                                        '-o', os.path.join(self.blobtools_dir, sample_name, sample_name),
                                        '--nodes', nodes,
                                        '--names', names
                                        ])
            else:
                exitcode = run_command([
                                        "blobtools", 'create',
                                        '-i', assembly, 
                                        '-b', os.path.join(self.bam_dir, '{}.bam'.format(sample_name)),
                                        '-t', blast_hits_file,
                                        '-x', 'bestsumorder',
                                        '-o', os.path.join(self.blobtools_dir, sample_name, sample_name),
                                        '--nodes', nodes,
                                        '--names', names
                                        ])
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'blobtools create', sample_name)
    
    def run_bt_view(self, reads):
        sample_name = get_sample_name(reads)
        blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
        blob_table = os.path.join(self.blobtools_dir, sample_name, '{}.{}.blobDB.table.txt'.format(sample_name,sample_name))
        if not os.path.isfile(blob_table):
            logging.info("{}: blobtools: Generating blob table...".format(sample_name))
            exitcode = run_command([
                                    "blobtools", 'view',
                                    '-i', blobdb_file,
                                    '-x', 'bestsumorder',
                                    '-r', 'species',
                                    '-b',
                                    '-o', os.path.join(self.blobtools_dir, sample_name, sample_name)
                                    ])
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'blobtools view', sample_name)
    
    def run_bt_plot_genus(self, reads):
        sample_name = get_sample_name(reads)
        blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
        blob_plot_genus = os.path.join(self.blobtools_dir, sample_name, '{}.{}.blobDB.json.bestsumorder.genus.p7.span.100.blobplot.spades.png'.format(sample_name, sample_name))
        if not os.path.isfile(blob_plot_genus):
            logging.info("{}: blobtools: Generating blob genus plot...".format(sample_name))
            exitcode = run_command([
                                    "blobtools", 'plot',
                                    '-i', blobdb_file,
                                    '-x', 'bestsumorder',
                                    '-r', 'genus',
                                    '-l', '100',
                                    '-o', os.path.join(self.blobtools_dir, sample_name, sample_name)
                                    ])
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'blobtools plot genus', sample_name)
    
    def run_bt_plot_species(self, reads):
        sample_name = get_sample_name(reads)
        blobdb_file = os.path.join(self.blobtools_dir, sample_name, '{}.blobDB.json'.format(sample_name))
        blob_plot_species = os.path.join(self.blobtools_dir, sample_name, '{}.{}.blobDB.json.bestsumorder.species.p7.span.100.blobplot.spades.png'.format(sample_name, sample_name))
        if not os.path.isfile(blob_plot_species):
            logging.info("{}: blobtools: Generating blob species plot...".format(sample_name))
            exitcode = run_command([
                                    "blobtools", 'plot',
                                    '-i', blobdb_file,
                                    '-x', 'bestsumorder',
                                    '-r', 'species',
                                    '-l', '100',
                                    '-o', os.path.join(self.blobtools_dir, sample_name, sample_name)
                                    ])
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'blobtools plot species', sample_name)
    
    def run_prokka(self, reads):
        sample_name = get_sample_name(reads)
        make_outdirs(self.ann_ass_dir)
        assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
        output =  os.path.join(self.ann_ass_dir, sample_name, '{}.gbf'.format(sample_name))
        if not os.path.isfile(output):
            logging.info("{}: Annotating assembly with prokka...".format(sample_name))
            exitcode = run_command([
                                    'prokka',
                                    '--outdir', os.path.join(self.ann_ass_dir, sample_name),
                                    '--prefix', sample_name,
                                    '--cpus', self.threads,
                                    assembly
                                    ])
        else:
            exitcode = 0
        self.log_exit_status(exitcode, 'prokka', sample_name)
    
    def collate_sequence_files(self, reads):
        sample_name = get_sample_name(reads)
        make_outdirs(self.final_fasta)
        make_outdirs(self.final_gbk)
        make_outdirs(self.protein_fasta)
        fasta_assembly = os.path.join(self.raw_ass_dir, sample_name, 'scaffolds.fasta')
        gb_assembly = os.path.join(self.ann_ass_dir, sample_name, '{}.gbf'.format(sample_name))
        protein_fasta = os.path.join(self.ann_ass_dir, sample_name, '{}.faa'.format(sample_name))
        output_fna =  os.path.join(self.final_fasta, '{}.fna'.format(sample_name))
        output_gbk =  os.path.join(self.final_gbk, '{}.gb'.format(sample_name))
        output_faa =  os.path.join(self.protein_fasta, '{}.faa'.format(sample_name))
        if os.path.isfile(fasta_assembly) and not os.path.isfile(output_fna):
            shutil.copyfile(fasta_assembly, output_fna)
        if os.path.isfile(gb_assembly) and not os.path.isfile(output_gbk):
            shutil.copyfile(gb_assembly, output_gbk)
        if os.path.isfile(protein_fasta) and not os.path.isfile(output_faa):
            shutil.copyfile(protein_fasta, output_faa)
    
    
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
        self.run_fastp()
        self.assemble_reads_spades()
        self.run_quast()
        make_outdirs(self.bam_dir)
        for reads in self.reads_list:
            if get_sample_name(reads) in self.successful['assembly']:
                if self.searcher == 'blastn':
                    self.run_blastn(reads)
                if self.searcher == 'diamond':
                    self.run_diamond(reads)
                self.run_bwa_samtools(reads)
                self.run_bt_create(reads)
                self.run_bt_view(reads)
                self.run_bt_plot_genus(reads)
                self.run_bt_plot_species(reads)
        for reads in self.reads_list:
            if get_sample_name(reads) in self.successful['assembly']:
                self.run_prokka(reads)
                self.collate_sequence_files(reads)
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
        logging.debug(line.decode('utf-8').strip())


def make_outdirs(dir):
    if os.path.isdir(dir) == False:
        os.makedirs(dir)


def run_command(cmd):
    logging.info("Executing command: {}".format(" ".join(cmd)))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with process.stdout as std_out:
        log_subprocess_output(std_out)
    exitcode = process.wait()
    return exitcode


def run_parallel_commands(cmds, program, threads):
    threads = str(int(float(threads)))
    logging.info("Executing {} instances of {} in parallel".format(threads, program))
    process = subprocess.Popen(['parallel', '-j', threads], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    output, error = process.communicate(input='\n'.join(cmds))
    logging.debug(output)
    return process.returncode


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

def run(threads, mem, reads_dir, force, outdir, verbosity, db, ion_torrent, assembler, searcher):
    configure_outdir(outdir, force)
    configure_logs(verbosity, outdir)
    logging.info("Beginning asqcan run...")
    job = Asqcan(threads, mem, db, reads_dir, outdir, verbosity, ion_torrent, assembler, searcher)
    if force == True:
        logging.info("Force option given. Deleting previous asqcan run...")
    logging.info("Log for the current run can be found at {}.".format(os.path.join(outdir, 'logs')))
    job.run_asqcan()
    logging.info("Run complete. Thanks for using asqcan.")
    job.clear_temp_dir()

