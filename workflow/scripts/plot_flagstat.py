import sys
import re
from matplotlib import pyplot as plt

samples = []
mapped = []

for flagstat in sys.argv:
    if flagstat.endswith(".flagstat"):
        with open(flagstat) as file:
            for line in file:
                x = re.findall("mapped \([0-9]+[.]{1}[0-9]+[%]{1}", line)
                if len(x) != 0:
                    samples.append(flagstat.split('/')[-1].split('.')[0])
                    mapped.append(float(re.findall("[0-9]+[.]{1}[0-9]+", x[0])[0]))

plt.bar(samples, mapped)
plt.title("% of mapped reads for each sample")
plt.xlabel("Samples")
plt.ylabel("% of mapped reads")
plt.savefig(sys.argv[-1], format='png')
