epiMuller README
================

.. figure:: https://raw.githubusercontent.com/jennifer-bio/epiMuller/main/images/case_scaled_lineages_long.png
   :alt: Muller plot image

   Muller plot image

About
-----

Author
~~~~~~

Jennifer L Havens

Purpose
~~~~~~~

Visualize lineages overtime, with phylogentic context, based on viral
genomes

Language
~~~~~~~~

Python3

Inputs
~~~~~~

Alingment, collection date, PANGO lineage, Nextstain JSON files, and
timetree

Source code avaliblity
~~~~~~~~~~~~~~~~~~~~~~

`gitHub <https://github.com/jennifer-bio/epimuller>`__

Documentation avaliblity
~~~~~~~~~~~~~~~~~~~~~~~~

`Read the Docs <https://epimuller.readthedocs.io/en/stable/>`__

Quick start
-----------

::

    pip3 install epimuller

    epimuller [-h] [-oDir OUTDIRECTORY] -oP OUTPREFIX -n
         INNEXTSTRAIN -m INMETA [-p INPANGOLIN]
         [-f TRAITOFINTERSTFILE]
         [-k TRAITOFINTERSTKEY]
         [-aa AAVOCLIST [AAVOCLIST ...]]
         [-t TIMEWINDOW] [-s STARTDATE] [-e ENDDATE]
         [-mt MINTIME] [-min MINTOTALCOUNT]
         [-c CASES_NAME] [-l {date,time}]
         [-lp {Right,Max,Start,End}]

SOME EXAMPLES
-------------

Examples for full run
~~~~~~~~~~~~~~~~~~~~~

To prep files for these examples for epimuller look at
Example\_CommandsFromScratch.txt

Visulize default aa mutation list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    epimuller \
        -n inputData/GISAID_NYCPHL_04_29/02_nextstrainResults \
        -m inputData/GISAID_NYCPHL_04_29/gisaid_2021_04_30_00_rename.tsv \
        -oDir 03_results_NYCPHL_April29 \
        -oP 01_defaultAAList \
        -c inputData/CITY_US-NY_NYC_outbreakinfo_epidemiology_data_2021-04-30.tsv

Visulize a trait: lineage
^^^^^^^^^^^^^^^^^^^^^^^^^

::

    epimuller \
        -n inputData/GISAID_NYCPHL_04_29/02_nextstrainResults \
        -m inputData/GISAID_NYCPHL_04_29/gisaid_2021_04_30_00_rename.tsv \
        -oDir 03_results_NYCPHL_April29 \
        -oP 02_pangolin \
        -c inputData/CITY_US-NY_NYC_outbreakinfo_epidemiology_data_2021-04-30.tsv \
        --traitOfInterstFile traits.json \
        --traitOfInterstKey lineage \
        -lp Max \
        -min 100 \

Visulize your own aa mutation list
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    epimuller \
        -n inputData/GISAID_NYCPHL_04_29/02_nextstrainResults \
        -m inputData/GISAID_NYCPHL_04_29/gisaid_2021_04_30_00_rename.tsv \
        -oDir 03_results_NYCPHL_April29 \
        -oP 03_selectedAA \
        -c inputData/CITY_US-NY_NYC_outbreakinfo_epidemiology_data_2021-04-30.tsv \
        -aa 'SE484K' 'S*452*' \
        -min 50 \ 
        -mt 20

Known edge cases / featrues to add
----------------------------------

Known edge cases which are not correctly dealt with or features I intend
to add (that I will get around to fixing eventually) If you run into
anything else please let me know on
`gitHub <https://github.com/jennifer-bio/epimuller>`__

::

        - nt_muts ; not set up for nt mutations (only amino acid or trait)
        - only takes nextstrain json files - intending to set up to take treetime output
        - feel free to ignore the undefined.svg that gets made - it is related to checking the size of the text to space out labels
        - add plot and font size to arg options

Addtional features
------------------

Color
~~~~~

If you would like to specify color for clade: in --parentHierarchy\_name
file (of epimuller-draw/drawMuller.py input) add col with name: "color"
and hex color value (starting with #) for clades you want to specify.

Plot and font size
~~~~~~~~~~~~~~~~~~

In the file: scripts/drawMuller.py ; near top of script change value for
desired WIDTH, HEIGHT, LEGENDWIDTH (space on right side of plot for
labels), MARGIN, or FONTSIZE variables. Then run with source code in
epimuller folder by $python3 scripts/drawMuller.py [Arugments]

Parse GISAID fasta for metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

epimuller-parse If you have downloaded sequences from GISAID under the
search tab, you can parse out the names into a metadata file (format
tested as of 2021-04-30)

ARGUMENTS
---------

::

    optional arguments:
        -h, --help            show this help message and exit

    Options for full repot:
        -oDir OUTDIRECTORY, --outDirectory OUTDIRECTORY
            folder for output (default: ./)
        -oP OUTPREFIX, --outPrefix OUTPREFIX
            prefix of out files withen outDirectory (default:
            None)

    Options passed to epimuller-define:
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

    Options passed to epimuller-draw:
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

Only make abundance and hiearchy files
--------------------------------------

::

    usage: epimuller-define [-h] -n INNEXTSTRAIN -m INMETA [-p INPANGOLIN]
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

Only plot
---------

::

    usage: epimuller-draw [-h] -p PARENTHIERARCHY_NAME -a ABUNDANCE_NAME
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

Install methods
---------------

With pip
~~~~~~~~

::

    pip3 install epimuller

    #If there is an issue with cairo, try:

    pip3 install pycairo
    pip3 install epimuller

From source
~~~~~~~~~~~

Download source code from
`gitHub <https://github.com/jennifer-bio/epimuller>`__ or
`pypi <https://pypi.org/project/epimuller/>`__

::

    #open as needed for download format
    tar -zxvf epimuller-[version].tar.gz

    cd epimuller-[version]

    python3 setup.py install

With Bioconda
~~~~~~~~~~~~~

Does not work at time of writing

::

    #conda install -c bioconda epimuller

Run scripts directly
~~~~~~~~~~~~~~~~~~~~

This is currently the best way to change display and font size, howerver
you will have to install all dependencies.

Download source code from
`gitHub <https://github.com/jennifer-bio/epimuller>`__ or
`pypi <https://pypi.org/project/epimuller/>`__

::

    #open as needed for download format
    tar -zxvf epimuller-[version].tar.gz

    cd epimuller-[version]

    #to run epimuller
    python3 ./scripts/mutationLinages_report.py [arugments]

    #to run epimuller-parse
    python3 ./scripts/parseFastaNames.py  [arugments]

    #to run epimuller-define 
    python3 ./scripts/defineAndCountClades.py  [arugments]

    #to run epimuller-draw 
    python3 ./scripts/drawMuller.py  [arugments]

Citation
--------

Please `link to this
github <https://github.com/jennifer-bio/epimuller>`__ if you have used
epimuller in your research.

Extra notes on GISAID
~~~~~~~~~~~~~~~~~~~~~

If you do use GISAID data please acknowledge the contributers, such as
with `language suggested by
GISAID <https://www.gisaid.org/help/publish-with-data-from-gisaid/>`__.
