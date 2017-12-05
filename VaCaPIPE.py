from chunkypipes.components import *

class Pipeline(BasePipeline):
	def dependencies(self):
		return ['cutadapt', 'pandas']

	def description(self):
		return 'Pair-end sequencing variant call pipeline for reads > 70bp'

	def configure(self):
		return {
			'annovar':{
				'path': 'Full path to ANNOVAR perl scripts ',
				'genome_build':'Genome build to use for annotation (see annovar documentation for options)'
			},
			'bwa':{
				'path': 'Full path to BWA executable',
				'threads': 'Number of threads to use for bwa mem',
				'reference': 'Directory reference genome in FASTA format'
			},
			'cutadapt':{
				'path': 'Full path to cutadapt executable',
				'fwd_three_prime': 'Adapater sequences for 3\' forward strand (leave blank if NONE)',
				'rev_three_prime': 'Adapater sequences for 3\' reverse strand (leave blank if NONE)',
				'fwd_five_prime': 'Adapater sequences for 5\' forward strand (leave blank if NONE)',
				'rev_five_prime': 'Adapater sequences for 5\' reverse strand (leave blank if NONE)',
				'base_quality_score_min': 'Minimum quality score to be trimmed at both the 3\' and 5\' end of reads - assumes Phred +33 scoring',
				'read_length_min': 'Minimum read length before both pairs are discarded'
			},
			'databases':{
				'path': 'Full path to directory of all database files to use for ANNOVAR'
			}
			'fastqc':{
				'path':'Full path to FastQC executable'
			},
			'platypus':{
				'path': 'Full path to Playtpus executable'
			},
			'samtools':{
				'path':'Full path to samtools executable'
			},
			'vcftools':{
				'path': 'Full path to vcftools executable'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('--fastq', type=str, help='tab-delimited file with list of FASTQ files to process (full path).  Pairs need to be on same line separated by tab.  One pair per line.  Files should have .fastq.gz or .fq.gz extension')
		parser.add_argument('--outDir', default=os.getcwd(), type=str, help='Full directory path to output files')
		parser.add_argument('--startStep', default='rawSeq_qc', type=str)
		parser.add_argument('--endStep', type=str)
		parser.add_argument('--parallelizeJobs', default=1, type=int, help='number of jobs to run concurrently/parallelize across n CPUs.  Must be integer.  Default:1')
		parser.add_argument('--annotationDB', default=None, type=str, help='comma-separated list of databases to use for ANNOVAR; should be same format as -protocol function in annovar table_annovar.pl')
		parser.add_argument('--annotationType', default=None, type=str, help='omma-separated list of database type to use for ANNOVAR; should be same format as -operation function in annovar table_annovar.pl')

	def run_pipeline(self, pipeline_args, pipeline_config):
		import os
		import pandas
		import subprocess
		import re


		@staticmethod
		def check_steps():
			['rawSeq_qc', 'alignment', 'variantCall', 'filtering', 'annotation']
			pass

		# check adapters that are blank
		if pipeline_config['cutadapt']['fwd_three_prime'] == '':
			pipeline_config['cutadapt']['fwd_three_prime'] = 'ZZZZZZ'
		if pipeline_config['cutadapt']['rev_three_prime'] == '':
			pipeline_config['cutadapt']['rev_three_prime'] = 'ZZZZZZ'
		if pipeline_config['cutadapt']['fwd_five_prime'] == '':
			pipeline_config['cutadapt']['fwd_five_prime'] = 'ZZZZZZ'
		if pipeline_config['cutadapt']['rev_five_prime'] == '':
			pipeline_config['cutadapt']['rev_five_prime'] = 'ZZZZZZ'

		
		# create shortcut to executable file
		fastqc = Software('fastqc', pipeline_config['fastqc']['path'])
		cutadapt = Software('cutadapt', pipeline_config['cutadapt']['path'])
		bwa_mem = Software('bwa', pipeline_config['bwa']['path'] + ' mem')
		bwa_idx = Software('bwa', pipeline_config['bwa']['path'] + ' index')
		flagstats = Software('samtools', pipeline_config['samtools']['path'] +  ' flagstats')
		view = Software('samtools', pipeline_config['samtools']['path'] + ' view')
		sort = Software('samtools', pipeline_config['samtools']['path'] +  ' sort')
		index = Software('samtools', pipeline_config['samtools']['path'] + ' index')

		
		if step[0] == 'rawSeq_qc':
			get_loc = re.compile('^(.*)(\.fastq\.gz|\.fq\.gz)')
			with open(pipeline_args['fastq']) as fastq_files:
				with ParallelBlock(processes = pipeline_args['parallelizeJobs']) as pblock:
					trimmed_files = []
					for fastq in fastq_files:
						fastq = fastq.strip().split('\t')
						fastq_fwd_name = get_loc.match(fastq[0])
						fastq_rev_name = get_loc.match(fastq[1])
						trimmed_files.append((fastq_fwd_name + '_trimmed.fastq.gz', fastq_rev_name + '_trimmed.fastq.gz'))
						pblock.add(
							cutadapt.prep(
								Parameter('-a', pipeline_config['cutadapt']['fwd_three_prime']),
								Parameter('-A', pipeline_config['cutadapt']['rev_three_prime']),
								Parameter('-g', pipeline_config['cutadapt']['fwd_five_prime']),
								Parameter('-G', pipeline_config['cutadapt']['rev_five_prime']),
								Parameter('-q', pipeline_config['cutadapt']['base_quality_score_min'] + ',' + pipeline_config['cutadapt']['base_quality_score_min']),
								Parameter('--pair-filter=any'),
								Parameter(fastq_fwd_name.group(1) + '_trimmed.fastq.gz'),
								Parameter('-p', fastq_rev_name.group(1) + '_trimmed.fastq.gz'),
								Parameter(fastq[0]),
								Parameter(fastq[1])
							)
						)
						
						pblock.add(
							fastqc.prep(
								Parameter(fastq_fwd_name.group(1) + '_trimmed.fastq.gz'),
								Parameter(fastq_rev_name.group(1) + '_trimmed.fastq.gz')
							)
						)

			step.pop(0)

			
		if step[0] == 'alignment':

			# check if reference genome has bwa index file
			get_loc = re.compile('^(.*[\\\/])/(.*)') # get directory where ref is located -> capture group 1; name of ref_file.fa -> capture group 2
			loc = get_loc.match(pipeline_config['bwa']['reference'])
			if loc.group(1) + '/' loc.group(2) + '.idx' in os.listdir(loc.group(1)):
				continue;
			else:
				bwa_idx.run(
					Parameter(pipeline_config['bwa']['reference'])
				)


			if pipeline_args['startStep'] != 'rawSeq_qc':
				trimmed_files = []
				with open() as trimmed_files:
					pass
					#TODO

			
			with ParallelBlock(processes = pipeline_args['parallelizeJobs']) as pblock:
				var_call_files = []
				get_prefix = re.compile('(.*)(\.fastq\.gz|\.fq\.gz)') # gets name of fastq file name prefix -> capture group 1; fastq.gz or fq.gz extension -> capture group 2
				for pairs in trimmed_files:
					fwd = get_prefix.match(pairs[0])
					rev = get_prefix.match(pairs[1])
					pblock.add(
						bwa_mem.prep(
							Parameter(pipeline_config['bwa']['reference']),
							Parameter('-t', pipeline_config['bwa']['threads']),
							Parameter(pairs[0]),
							Parameter(pairs[1]),
							Pipe( # this pipe will sort and create BAM bypass SAM file
								sort.prep(
									Parameter('-O', 'BAM'),
									Parameter('-o', str(fwd[1])+'_'+rev[1]+'.bam'),
									Parameter('-'),	
								))
					))

					# index bam file
					index.run(
						Parameter(str(fwd[1])+'_'+str(rev[1])+'.bam')
						)
					
					# record flagestats on bamfile
					flagstats.run(
						Parameter(str(fwd[1])+'_'+str(rev[1])+'.bam'),
						Redirect(stream=Redirect.STDOUT, dest=str(fwd[1])+'_'+str(rev[1])+'.flagstats')
						)

					var_call_files.append(str(fwd[1])+'_'+rev[1]+'.bam')
 
			step.pop(0)

		if step[0] == 'variantCall':

			step.pop(0)

		if step[0] == 'filtering':
			step.pop(0)


		if step[0] == 'annotation':
			step.pop[0]


