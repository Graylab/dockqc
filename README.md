# BCAPIN #

![BCAPIN](./bcapin.jpg)

Official repository for Benchmark of CArbohydrate Protein INteractions:

# THIS GITHUB REPO IS IN ALPHA. ISSUES AND QOL MAY BE MISSING AT THIS TIME. #
## WE PLAN TO RESOLVE ALL ISSUES AND QOL IMPLEMENTATION ISSUES In SEPTEMBER 2025 ##

```
Code writer (SWC) defended his dissertation and is on vacation and will return and make the code more up to his standard upon his return. We apologize for any inconvenience and please contact scanner1 [@] jhu.edu if you have any difficulty with the repo and its usage
```

Current issues:
```
.cif files integration
improved documentation
```

## Predicted structure Data

All generated structural data is available at `./bcapin_all_data.zip`

## DockQC Installation

### Local Install
```
conda env create -f bcapin.yml
conda activate bcapin
```

# How to run: Command Line #


```
python dockqc.py -x path/to/exp_file.pdb -p path/to/pred_file.pdb -xp A -xc B -pp A -pc B
```

### Help info ###
```
-x experimentally solved structure
-p predicted get_structure
-xp protein chain in experimentally solved structure (comma seperated)
-xc carbohydrate chain in experimentally solved structure (comma seperated)
-pp protein chain in predicted structure (comma seperated)
-pc carbohydrate chain in predicted structure (comma seperated)
```

### Example input ###
For `7EQR` example case provided here

```
python dockqc.py -x ex/wt.pdb -p ./ex/pred.pdb -xp B -xc H -pp A -pc B
```
