Using MicroPITA commandline

These common commands can be used on the default data set obtained when downloading MicroPITA, simply cut and paste them into a commandline in the downloaded micropita directory.


A. Expected input file.

I. PCL file
Although some defaults can be changed, microPITA expects a PCL file as an input file. A PCL file is a TEXT delimited file similar to an excel spread sheet with the following characteristics.

1. Rows represent metadata and features (bugs), columns represent samples.
2. The first row by default should be the sample ids.
3. Metadata rows should be next.
4. Lastly, rows containing features (bugs) measurements (like abundance) should be after metadata rows.
5. The first column should contain the ID describing the column. For metadata this may be, for example, "Age" for a row containing the age of the patients donating the samples. For measurements, this should be the feature name (bug name).
5. By default the file is expected to be TAB delimited.
6. If a consensus lineage or hierachy of taxonomy is contained in the feature name, the default delimiter between clades is the pipe ("|").

II. Targeted feature file
If using the targeted feature methodology, you will need to provide a txt file listing the feature(s) of interest. Each feature should be on it's own line and should be written as found in the input PCL file.


B. Basic unsupervised methods.

There are four unsupervised methods which can be performed:
diverse (maximum diversity), extreme (most dissimilar), representative (representative dissimilarity) and features (targeted feature).

The first three methods are performed as follows (selecting a default 10 samples):

$ python MicroPITA.py --lastmeta Label -m representative input/Test.pcl output.txt

$ python MicroPITA.py --lastmeta Label -m diverse input/Test.pcl output.txt

$ python MicroPITA.py --lastmeta Label -m extreme input/Test.pcl output.txt

Each of the previous methods are made up of the following pieces:
1. python MicroPITA.py to call the MicroPITA script.
2. --lastmeta which indicates the keyword (first column value) of the last row that contains metadata.
3. -m which indicates the method to use in selection.
4. input/Test.pcl which is the first positional argument indicating an input file
5. output.txt which is the second positional argument indicating the location to write to the output file.

Selecting specific features has additional arguments to consider --targets (required) and --feature_method (optional).

$ python MicroPITA.py --lastmeta Label -m features --targets input/TestFeatures.taxa input/Test.pcl output.txt

$ python MicroPITA.py --lastmeta Label -m features --feature_method abundance --targets input/TestFeatures.taxa input/Test.pcl output.txt

These additional arguments are described as:
1. --targets The path to the file that has the features (bugs or clades) of interest. Make sure they are written as they appear in your input file!
2. --feature_method is the method of selection used and can be based on ranked abundance ("rank") or abundance ("abundance"). The default value is rank.
To differentiate the methods, rank tends to select samples in which the feature dominates the samples reguardless of it's abundance.
Abundance tends to select samples in which the feature is most abundant without a guarentee that the feature is the most abundant feature in the sample. 


C. Basic supervised methods.

Two supervised methods are also available:
distinct and discriminant

These methods require an additional argument --label which is the first column keyword of the row used to classify samples for the supervised methods.
These methods can be performed as follows:

$ python MicroPITA.py --lastmeta Label --label Label -m distinct input/Test.pcl output.txt

$ python MicroPITA.py --lastmeta Label --label Label -m discriminant input/Test.pcl output.txt


D. Changing defaults.

Sample Selection:
To change the number of selected samples for any method use the -n argument. This example selects 6 representative samples instead of the default 10.

$ python MicroPITA.py --lastmeta Label -m representative -n 6 input/Test.pcl output.txt

When using a supervised method this indicates how many samples will be selected per class of sample. For example if you are performing supervised selection of 6 samples (-n 6) on a dataset with 2 classes (values) in it's label row, you will get 6 x 2 = 12 samples. If a class does not have 6 samples in it, you will get the max possible for that class. In a scenario where you are selecting 6 samples (-n 6) and have two classes but one class has only 3 samples then you will get 6 + 3 = 9 selected samples.

Stratification:
To stratify any method use the --stratify argument which is the first column keyword of the metadata row used to stratify samples before selection occurs. (Selection will occur independently within each strata). This example stratifies diverse selection by the "Label".

$ python MicroPITA.py --lastmeta Label --stratify Label -m representative input/Test.pcl output.txt

$ python MicroPITA.py --lastmeta Label --label Label --stratify StratifyLabel -m distinct input/Test.pcl output.txt

Sample ID header:
MicroPITA assumes the first row of the input file is the sample IDs, if it is not you may use --id to indicate the row.
--id expects the entry in the first column of your input file that matches the row used as Sample Ids. See the input file and the following command as an example.

$ python MicroPITA.py --id Sample --lastmeta Label -m representative input/Test.pcl output.txt

MicroPITA assumes the input file is TAB delimited, we strongly recommend you use this convention. If not, you can use --delim to change the delimiter used to read in the file.
Here is an example of reading the comma delimited file micropita/input/CommaDelim.pcl

$ python MicroPITA.py --delim , --lastmeta Label -m representative input/CommaDelim.pcl output.txt

MicroPITA assumes the input file has feature names in which, if the name contains the consensus lineage or full taxonomic hierarchy, it is delimited with a pipe "|". We strongly recommend you use this default. The delimiter of the feature name can be changed using --featdelim. Here is an example of reading in a file with periods as the delimiter.

$ python MicroPITA.py --featdelim . --lastmeta Label -m representative input/PeriodDelim.pcl output.txt


This covers how to use microPITA. Thank you for using this software and good luck with all your endeavors!
