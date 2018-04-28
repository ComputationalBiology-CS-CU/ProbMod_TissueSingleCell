import numpy as np
import matplotlib.pyplot as plt
from collections import Counter






if __name__ == "__main__":
	################################################
	filename = "../out/out_10.vcf"
	file = open(filename, 'r')
	count = 0
	upper = 30							## NOTE: we have redundant lines
	while count<upper:
		file.readline()
		count+=1

	##
	print("+++++++++++++++++++++++ id line:")
	print(file.readline().strip())

	##
	print("+++++++++++++++++++++++ content:")
	repo = Counter()
	list_count = []
	count00 = 0
	count01 = 0
	count11 = 0
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		##
		line0 = line[::]
		line = line.split('\t')[9:]				## NOTE: the position to start
		line = list(map(lambda x: x.split(':')[0], line))

		## number of genotyped individuals
		count = len(line) - line.count("./.")
		if count == 1:
			continue
		list_count.append(count)

		##
		repo += Counter(line)

		## NOTE we treat all heterozygous as 0/1, and all alternative homozygous as 1/1
		count00 += line.count('0/0')
		count01 += line.count('0/1')
		count01 += line.count('0/2')
		count11 += line.count('1/1')
		count11 += line.count('1/2')
		count11 += line.count('2/2')
		count11 += line.count('2/3')
	file.close()
	print("genotype statistics (all sites all individuals):")
	print(repo)
	print("count00, count01, count11:")
	print(count00, count01, count11)

	## histogram
	list_count = np.array(list_count)
	n, bins, patches = plt.hist(list_count, bins=16, alpha=0.75)
	#plt.axvline(list_count.mean(), color='k', linestyle='dashed', linewidth=1)
	#
	plt.xlabel('# of non-./. genotype')
	plt.ylabel('frequency')
	plt.title('Hist for # of genotyped individuals (remove 1 count)')
	plt.grid(True)
	plt.show()



