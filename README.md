# epiMuller

![Muller plot](https://raw.githubusercontent.com/jennifer-bio/epiMuller/main/images/case_scaled_lineages_long.png)

##### Author: 
Jennifer L Havens
##### Purpose: 
Visualize frequency of SARS-CoV2 variants, and with phylogenetic context, based on sequencing data, over time using muller plot
##### Language: 
Python3
##### Inputs: 
Alingment, collection date, PANGO lineage, Nextstain JSON files, and timetree

## Quick start

Best way to run is to run python 00_scripts/mutationLinages_report.py from command line in muller_plotting folder:


cd muller_plotting

```
python 00_scripts/mutationLinages_report.py [-h] [-oDir OUTDIRECTORY] -oP OUTPREFIX -n
                                 INNEXTSTRAIN -m INMETA [-p INPANGOLIN]
                                 [-f TRAITOFINTERSTFILE]
                                 [-k TRAITOFINTERSTKEY]
                                 [-aa AAVOCLIST [AAVOCLIST ...]]
                                 [-t TIMEWINDOW] [-s STARTDATE] [-e ENDDATE]
                                 [-mt MINTIME] [-min MINTOTALCOUNT]
                                 [-c CASES_NAME] [-l {date,time}]
                                 [-lp {Right,Max,Start,End}]

```

## SOME EXAMPLES 

```
python 00_scripts/mutationLinages_report.py -n 01_inputData -m 01_inputData/NYC_PHL_PRL_2021_03_24_ref.tsv -oDir 02_results -oP 01_defaultAAList -min 50 -c 01_inputData/CITY_US-NY_NYC_outbreakinfo_epidemiology_data_2021-04-23.tsv


python 00_scripts/mutationLinages_report.py -n 01_inputData -m 01_inputData/NYC_PHL_PRL_2021_03_24_ref.tsv -oDir 02_results -oP 02_traits -c 01_inputData/CITY_US-NY_NYC_outbreakinfo_epidemiology_data_2021-04-23.tsv --traitOfInterstFile temp_subclades.json --traitOfInterstKey clade_membership -l time -mt 40


python 00_scripts/mutationLinages_report.py -n 01_inputData -m 01_inputData/NYC_PHL_PRL_2021_03_24_ref.tsv -oDir 02_results -oP 03_allS_452mut -c 01_inputData/CITY_US-NY_NYC_outbreakinfo_epidemiology_data_2021-04-23.tsv -aa 'S*452*' -min 5 -lp Start

python 00_scripts/mutationLinages_report.py  -n 01_inputData -m 01_inputData/NYC_PHL_PRL_2021_03_24_ref.tsv -oDir 02_results -oP 04_E484K -s 2021-01-01 -e 2021-03-20 -mt 1 -c 01_inputData/CITY_US-NY_NYC_outbreakinfo_epidemiology_data_2021-04-03.tsv -aa 'SE484K' -min 10 -lp Max
```



## Known edge cases / featrues to add  
Known edge cases which are not correctly dealt with or features I intend to add (that I will get around to fixing eventually) 
If you run into anything else please let me know (jhavens@ucsd.edu)

	- nt_muts ; not set up for nt mutations (only amino acid or trait)
	- only takes nextstrain json files - intending to set up to take treetime output
	- feel free to ignore the undefined.svg that gets made - it is related to checking the size of the text to space out labels

## Addtional features  
If you would like to specify color for clade: add col with name: "color" and hex color value (starting with #) for clades you want to specify.


## ARGUMENTS  

```
optional arguments:
  -h, --help            show this help message and exit

Options for full repot:
  -oDir OUTDIRECTORY, --outDirectory OUTDIRECTORY
                        folder for output (default: ./)
  -oP OUTPREFIX, --outPrefix OUTPREFIX
                        prefix of out files withen outDirectory (default:
                        None)

Options passed to defineAndCountClades.py:
  -n INNEXTSTRAIN, --inNextstrain INNEXTSTRAIN
                        nextstrain results with tree.nwk and
                        [traitOfInterst].json (default: None)
  -m INMETA, --inMeta INMETA
                        metadata tsv with 'strain' and 'date'cols, optional:
                        cols of trait of interst; and pangolin col named:
                        'lineage' or 'pangolin_lin' (default: None)
  -p INPANGOLIN, --inPangolin INPANGOLIN
                        pangolin output lineage_report.csv file, if argument
                        not supplied looks in inMeta for col with
                        'pangolin_lin' or 'lineage' (default: metadata)
  -f TRAITOFINTERSTFILE, --traitOfInterstFile TRAITOFINTERSTFILE
                        name of nextstrain [traitOfInterst].json in
                        'inNextstrain' folder (default: aa_muts.json)
  -k TRAITOFINTERSTKEY, --traitOfInterstKey TRAITOFINTERSTKEY
                        key for trait of interst in json file (default:
                        aa_muts)
  -aa AAVOCLIST [AAVOCLIST ...], --aaVOClist AAVOCLIST [AAVOCLIST ...]
                        list of aa of interest in form
                        [GENE][*ORAncAA][site][*ORtoAA] ex. S*501*, gaps
                        represed by X (default: None)
  -t TIMEWINDOW, --timeWindow TIMEWINDOW
                        number of days for sampling window (default: 7)
  -s STARTDATE, --startDate STARTDATE
                        start date in iso format YYYY-MM-DD or 'firstDate'
                        which sets start date to first date in metadata
                        (default: 2020-03-01)
  -e ENDDATE, --endDate ENDDATE
                        end date in iso format YYYY-MM-DD or 'lastDate' which
                        sets end date as last date in metadata (default:
                        lastDate)

Options passed to drawMuller.py:
  -mt MINTIME, --MINTIME MINTIME
                        minimum time point to start plotting (default: 30)
  -min MINTOTALCOUNT, --MINTOTALCOUNT MINTOTALCOUNT
                        minimum total count for group to be included (default:
                        10)
  -c CASES_NAME, --cases_name CASES_NAME
                        file with cases - formated with 'date' in ISO format
                        and 'confirmed_rolling' cases, in tsv format (default:
                        None)
  -l {date,time}, --xlabel {date,time}
                        Format of x axis label: ISO date format or timepoints
                        from start (default: date)
  -lp {Right,Max,Start,End}, --labelPosition {Right,Max,Start,End}
                        choose position of clade labels (default: Right)

```



## Only make abundance and hiearchy files 

```
usage: defineAndCountClades.py [-h] -n INNEXTSTRAIN -m INMETA [-p INPANGOLIN]
                               [-f TRAITOFINTERSTFILE] [-k TRAITOFINTERSTKEY]
                               [-aa AAVOCLIST [AAVOCLIST ...]]
                               [-oDir OUTDIRECTORY] -oP OUTPREFIX
                               [-t TIMEWINDOW] [-s STARTDATE] [-e ENDDATE]

optional arguments:
  -h, --help            show this help message and exit
  -n INNEXTSTRAIN, --inNextstrain INNEXTSTRAIN
                        nextstrain results with tree.nwk and
                        [traitOfInterst].json (default: None)
  -m INMETA, --inMeta INMETA
                        metadata tsv with 'strain' and 'date'cols, optional:
                        cols of trait of interst; and pangolin col named:
                        'lineage' or 'pangolin_lin' (default: None)
  -p INPANGOLIN, --inPangolin INPANGOLIN
                        pangolin output lineage_report.csv file, if argument
                        not supplied looks in inMeta for col with
                        'pangolin_lin' or 'lineage' (default: metadata)
  -f TRAITOFINTERSTFILE, --traitOfInterstFile TRAITOFINTERSTFILE
                        name of nextstrain [traitOfInterst].json in
                        'inNextstrain' folder (default: aa_muts.json)
  -k TRAITOFINTERSTKEY, --traitOfInterstKey TRAITOFINTERSTKEY
                        key for trait of interst in json file (default:
                        aa_muts)
  -aa AAVOCLIST [AAVOCLIST ...], --aaVOClist AAVOCLIST [AAVOCLIST ...]
                        list of aa of interest in form
                        [GENE][*ORAncAA][site][*ORtoAA] ex. S*501*, gaps
                        represed by X (default: None)
  -oDir OUTDIRECTORY, --outDirectory OUTDIRECTORY
                        folder for output (default: ./)
  -oP OUTPREFIX, --outPrefix OUTPREFIX
                        prefix of out files withen outDirectory (default:
                        None)
  -t TIMEWINDOW, --timeWindow TIMEWINDOW
                        number of days for sampling window (default: 7)
  -s STARTDATE, --startDate STARTDATE
                        start date in iso format YYYY-MM-DD or 'firstDate'
                        which is in metadata (default: 2020-03-01)
  -e ENDDATE, --endDate ENDDATE
                        end date in iso format YYYY-MM-DD or 'lastDate' which
                        is in metadata (default: lastDate)
```




## Only plot 

```
usage: drawMuller.py [-h] -p PARENTHIERARCHY_NAME -a ABUNDANCE_NAME
                     [-c CASES_NAME] -o OUTFOLDER [-mt MINTIME]
                     [-min MINTOTALCOUNT] [-l {date,time}]
                     [-lp {Right,Max,Start,End}]

optional arguments:
  -h, --help            show this help message and exit
  -p PARENTHIERARCHY_NAME, --parentHierarchy_name PARENTHIERARCHY_NAME
                        csv output from mutationLinages_report.py with child
                        parent col (default: None)
  -a ABUNDANCE_NAME, --abundance_name ABUNDANCE_NAME
                        csv output from mutationLinages_report.py with
                        abundances of clades (default: None)
  -c CASES_NAME, --cases_name CASES_NAME
                        file with cases - formated with 'date' in ISO format
                        and 'confirmed_rolling' cases, in tsv format (default:
                        None)
  -o OUTFOLDER, --outFolder OUTFOLDER
                        csv output from mutationLinages_report.py with child
                        parent col (default: None)
  -mt MINTIME, --MINTIME MINTIME
                        minimum time point to start plotting (default: 30)
  -min MINTOTALCOUNT, --MINTOTALCOUNT MINTOTALCOUNT
                        minimum total count for group to be included (default:
                        10)
  -l {date,time}, --xlabel {date,time}
                        Format of x axis label: ISO date format or timepoints
                        from start (default: date)
  -lp {Right,Max,Start,End}, --labelPosition {Right,Max,Start,End}
                        choose position of clade labels (default: Right)

```
