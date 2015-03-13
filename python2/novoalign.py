#!/usr/bin/env python

def main():

	#command line arguments
	parser = argparse.ArgumentParser(
		description = 'illumina qseq files > sam alignment file',
		epilog = '\a9 2011 Michael Yourshaw all rights reserved')
	parser.add_argument('--input', '-f', nargs='+',
											help='one or two paired-end qseq files files')
	parser.add_argument('--output', '-o',
											help='output sam file')
	parser.add_argument('--RG',
											help='read group ID')
	parser.add_argument('--CN',
											help='sequencing center (default: UCLA)')
	parser.add_argument('--DS',
											help='read group description')
	parser.add_argument('--DT',
											help='run date ISO8601M')
	parser.add_argument('--LB',
											help='library')
	parser.add_argument('--PI',
											help='predicted median insert size')
	parser.add_argument('--PL',
											help='platform (default=ILLUMINA)')
	parser.add_argument('--PU',
											help='platform unit (flowcell_lane[_barcode]')
	parser.add_argument('--SM',
											help='sample')
	parser.add_argument('--adapter1',
											help='')
	parser.add_argument('--adapter2',
											help='')
	parser.add_argument('--novoindex',
											help='')
	parser.add_argument('--config',
											help='configuration file (default=[pypeline.cfg,~/pypeline.cfg,pypeline.cfg])')
	parser.add_argument('--novoalign_options',
											help='novoalign options in the form "-opt value [...]"')
	parser.add_argument('--tmpdir',
											help='temporary directory (default from config->tmpscratch)')
	args = parser.parse_args()

	#Novoalign
	if [[ "${READGROUP}" == '' ]]; then rgid="${base%.pe1.qseq.txt}"; else rgid="${READGROUP}"; fi #GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA
	experiment="${rgid%%.*}"; #GMD108_XT_Illumina
	if [[ "${SAMPLE}" == '' ]]; then sample="${experiment%%_*}"; else sample="${SAMPLE}"; fi #GMD108
	if [[ "${LIBRARY}" == '' ]]; then library="${experiment%_*}"; else library="${LIBRARY}"; #GMD108_XT
	x="${rgid#*.}"; #2011-05-11.HWI-ST430.243.6.ACTTGAA
	if [[ "${DATE}" == '' ]]; then dt="${x%%.*}"; else dt="${DATE}"; #2011-05-11
	if [[ "${PLATFORMUNIT}" == '' ]]; then pu="${x#*.}"; else pu="${PLATFORMUNIT}"; #HWI-ST430.243.6.ACTTGAA

	rg='@RG\tID:'${rgid}'\tCN:UCLA\tDT:'${dt}'\tLB:'${library}'\tPL:ILLUMINA\tPU:'${pu}'\tSM:'${sample};
	#@RG\tID:GMD108_XT_Illumina.2011-05-11.HWI-ST430.243.6.ACTTGAA\tCN:UCLA\tDT:2011-05-11\tLB:GMD108_XT\tPL:ILLUMINA\tPU:HWI-ST430.243.6.ACTTGAA\tSM:GMD108
	
	name=novoalign_`basename ${sam}`;
	cmd="${NOVOALIGN} -k -o SAM \"$rg\" -d ${NOVOINDEX} -a ${ADAPTER1} ${ADAPTER2} -F QSEQ -f ${f} ${f2} > ${sam};";
	if [[ ${HIMEM} == 0 ]]; then NovoalignId=`echo "${cmd}" | ${qsub_lomem} -N ${name}` else NovoalignId=`echo "${cmd}" | ${qsub_himem} -N ${name}`; fi
	echo "${NovoalignId} ${name} $[cmd}";
	
if __name__ == "__main__": sys.exit(main())