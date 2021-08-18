These files contain the time series and the associated hashtags we obtained by sampling Twitter for our paper "Universality, 
criticality and complexity of information propagation on social media". The analysis is reported in 

LINK TO PAPER

Please acknowledge the use of these data by citing the paper above.


#################################
#################################


DATA ORGANIZATION


We created a single zip file with all the time series and a single zip file with all the hashtags. There is a one-to-one correspondence 
between lines in the two files. However, it was not possible to upload these two files as they are too heavy. We split the files and uploaded
numerous small files. The split was done via command line by

zip single_file.zip --out multiple_files.zip -s 25m

The data can be unzipped by first collecting them together and then unizpping the resulting file. This can be via command line by

zip -F multiple_files.zip --out single_file.zip
unzip single_file.zip

The first of these two lines recollect the files in a single one and the second unzip the resulting file. 


#################################
#################################


FILES CONTENT


As stated, here is a one-to-one correspondence between lines in the time series file and lines in the hashtags file, i.e., the hashtag 
stored in line X is the hashtag of the time series stored in line X. Time series are stored as follows:

Ka t1 t2 t3 \n
Kb t1 t2 t3 t4 t5 \n
.
.
.
Kn t1 t2 \n

where:

Ka, Kb,..., Kn is an integer specifying the number of events that compose the time series a, b,..., n respectively. In the example above 
we would have Ka=3, Kb=5, Kn=2.

t1 t2 ... is the time series, i.e., a sequence of chronologically ordered interevent times. The last interevent time, 
in our implementation, represents the distance between the end of the temporal window and the last event time. 
It thus does not represent an event. As stated in the Supplemental Material of our paper, the temporal window ranges from 
2019, October 1st to 2019, November 30th.
