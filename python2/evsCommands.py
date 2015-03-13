#!/bin/sh
#create shell script to download all data from http://snp.gs.washington.edu/EVS/
chrs = {
1:	249250621,
2:	243199373,
3:	198022430,
4:	191154276,
5:	180915260,
6:	171115067,
7:	159138663,
8:	146364022,
9:	141213431,
10:	135534747,
11:	135006516,
12:	133851895,
13:	115169878,
14:	107349540,
15:	102531392,
16:	90354753,
17:	81195210,
18:	78077248,
19:	59128983,
20:	63025520,
21:	48129895,
22:	51304566,
'X':	155270560,
'Y':	59373566
}
format = 'text' #'vcf'
for chr in sorted(chrs.keys()):
	start = 1
	last = chrs[chr]
	while start <= last:
		end = min(start+999999, last)
		region = '{}:{}-{}'.format(chr,start,end)
		print 'echo {}'.format(region)
		print 'java -jar /Users/myourshaw/apps/evsClient-v.0.0.1/evsClient.jar -t {} -f {} '.format(region, format)
		start = end +1
