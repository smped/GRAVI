import yaml
import pandas as pd
import sys

f = "config/config.yml"
config = yaml.safe_load(open(f))
df = pd.read_table(config['samples']['file'])

## Check for required column names
print("Checking " + config['samples']['file'])
req_cols = ['sample', 'target', 'treat', 'input']
if not set(req_cols).issubset(df.columns):
  print("Columns misspecified in " + config['samples']['file'])
  print("Can only be one of the following:\n")
  print(*req_cols, sep = ", ")
  sys.exit(1)

print(config['samples']['file'] + " has the correct column names\n")

## Now set all values as required
samples = list(set(df['sample']))
targets = list(set(df['target']))
treats = list(set(df['treat']))

#######################
## Check comparisons ##
#######################

print("Checking requested comparisons...")

## Ignore any targets requested which are not present in the data
comps = config['samples']['comparisons']
target_comps = set(comps.keys())
target_comps = target_comps.intersection(targets)

## Check all comparisons are valid
def check_comparisons(target):
  x = comps[target]
  ind = df['target'] == target
  if type(x[0]) is list:
    s = set()
    for i in range(len(x)):
      s.add(set(x[i]).issubset(df['treat'][ind]))
    return(all(s))
  if type(x[0]) is str:
    return(set(x).issubset(df['treat'][ind]))
    

chk_comps = [False] * len(target_comps)
for i in range(len(target_comps)):
  chk_comps[i] = check_comparisons(list(target_comps)[i])

if not all(chk_comps):
  for i in range(len(chk_comps)):
    if not(chk_comps[i]):
      print(
      "Treatment level mismatch for " + list(target_comps)[i] +
      "\nPlease check the requested treatments with those in " +
      config['samples']['file']
      )
  sys.exit(1)

print("Comparisons are specified correctly\n")
print("Checking covariates for IHW...")

## Check requested covariates for IHW. These can only be another TF or 
## logCPM/width
ihw = config['samples']['ihw']
target_ihw = set(ihw.keys())
target_ihw = target_ihw.intersection(targets)

for x in list(target_ihw):
  if not type(x) is str:
    print("IHW covariates not specified correctly for" + x)
    sys.exit(1)
  

for x in list(target_ihw):
  ind = df['target'] != x
  poss_vals = df['target'][ind]
  other_targets = set(ihw[x]).issubset(poss_vals)
  other_vals = set(ihw[x]).issubset(['logCPM', 'width'])
  if not ((other_targets and not other_vals) or other_vals):
    print("IHW covariates mis-specified for " + x)
    print(" Covariates can only be other targets or 'logCPM' or 'width'")
    sys.exit(1)
  else:
    print("IHW covariates appear specified correctly for " + x)

True
