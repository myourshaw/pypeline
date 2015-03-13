#!/usr/bin/env python

import copy

peds = {
	"mms2267":"mms 2267 0 0 1 1",
	"mms2268":"mms 2268 0 0 2 1",
	"mms2269":"mms 2269 2267 2268 1 2",
	"mms2270":"mms 2270 2267 2268 1 2",
	"mms2271":"mms 2271 2267 2268 1 2",
	"mms2272":"mms 2272 2267 2268 2 1",
	"mms2273":"mms 2273 2267 2268 1 2",
	"mms2274":"mms 2274 2267 2268 2 1",
	"mms2276":"mms 2276 2267 2268 1 1",
	"mms2277":"mms 2276 2267 2268 2 1",
	"mms2278":"mms 2276 2267 2268 2 1",
	"mms2488":"mms 2488 2267 2268 2 1"
}
pedcheck_peds = {
	"mms2267":"mms 2267 0 0 1",
	"mms2268":"mms 2268 0 0 2",
	"mms2269":"mms 2269 2267 2268 1",
	"mms2270":"mms 2270 2267 2268 1",
	"mms2271":"mms 2271 2267 2268 1",
	"mms2272":"mms 2272 2267 2268 2",
	"mms2273":"mms 2273 2267 2268 1",
	"mms2274":"mms 2274 2267 2268 2",
	"mms2276":"mms 2276 2267 2268 1",
	"mms2277":"mms 2276 2267 2268 2",
	"mms2278":"mms 2276 2267 2268 2",
	"mms2488":"mms 2488 2267 2268 2"
}
tr = {
	"AA":"1 1",
	"AB":"1 2",
	"BA":"1 2",
	"BB":"2 2",
	"NoCall":"0 0"
}
pedcheck_tr = {
	"AA":"1 1",
	"AB":"1 2",
	"BA":"1 2",
	"BB":"2 2",
	"NoCall":"0 0"
}
gt = []
gt_autosome = []
gt_x = []
out = open("/Users/myourshaw/lab/MMS/mathematica.mms/mms.ped","w")
out_autosome = open("/Users/myourshaw/lab/MMS/mathematica.mms/mms_autosome.ped","w")
out_x = open("/Users/myourshaw/lab/MMS/mathematica.mms/mms_x.ped","w")
#pedcheck_out = open("/Users/myourshaw/lab/MMS/mathematica.mms/pedcheck_mms.ped","w")
#pedcheck_out_autosome = open("/Users/myourshaw/lab/MMS/mathematica.mms/pedcheck_mms_autosome.ped","w")
#pedcheck_out_x = open("/Users/myourshaw/lab/MMS/mathematica.mms/pedcheck_mms_x.ped","w")
#SNP ID	Chromosome	Physical Position	dbSNP RS ID	mms2267	mms2267_conf	mms2268	mms2268_conf	mms2269	mms2269_conf	mms2270	mms2270_conf	mms2271	mms2271_conf	mms2272	mms2272_conf	mms2273	mms2273_conf	mms2274	mms2274_conf	mms2276	mms2276_conf	mms2488	mms2488_conf
for line in open("/Users/myourshaw/lab/MMS/mathematica.mms/mms.genotypes.decode.txt"):
	line = line.strip()
	if len(line) == 0:
		continue
	fields = line.split("\t")
	if fields[0] == "SNP ID":
		ids = fields[4:len(fields):2]
		continue
	gt += fields[4:len(fields):2]
	if fields[1].upper() == "X":
		gt_x += fields[4:len(fields):2]
	else:
		gt_autosome += fields[4:len(fields):2]
#pedcheck_gt = copy.deepcopy(gt)
#pedcheck_gt_autosome = copy.deepcopy(gt_autosome)
#pedcheck_gt_x = copy.deepcopy(gt_x)
for i in range(len(gt)):
	gt[i] = tr[gt[i]]
for i in range(len(gt_autosome)):
	gt_autosome[i] = tr[gt_autosome[i]]
for i in range(len(gt_x)):
	gt_x[i] = tr[gt_x[i]]
#for i in range(len(pedcheck_gt)):
#	pedcheck_gt[i] = pedcheck_tr[pedcheck_gt[i]]
#for i in range(len(pedcheck_gt_autosome)):
#	pedcheck_gt_autosome[i] = pedcheck_tr[pedcheck_gt_autosome[i]]
#for i in range(len(pedcheck_gt_x)):
#	pedcheck_gt_x[i] = pedcheck_tr[pedcheck_gt_x[i]]
for i in range(len(ids)):
	ped = peds[ids[i]]
	genotype = " ".join(gt[i:len(gt):len(ids)])
	out.write("%s %s\n" % (ped,genotype))
	genotype = " ".join(gt_autosome[i:len(gt_autosome):len(ids)])
	out_autosome.write("%s %s\n" % (ped,genotype))
	genotype = " ".join(gt_x[i:len(gt_x):len(ids)])
	out_x.write("%s %s\n" % (ped,genotype))
	#ped = pedcheck_peds[ids[i]]
	#genotype = " ".join(pedcheck_gt[i:len(pedcheck_gt):len(ids)])
	#pedcheck_out.write("%s %s\n" % (ped,genotype))
	#genotype = " ".join(pedcheck_gt_autosome[i:len(pedcheck_gt_autosome):len(ids)])
	#pedcheck_out_autosome.write("%s %s\n" % (ped,genotype))
	#genotype = " ".join(pedcheck_gt_x[i:len(pedcheck_gt_x):len(ids)])
	#pedcheck_out_x.write("%s %s\n" % (ped,genotype))
