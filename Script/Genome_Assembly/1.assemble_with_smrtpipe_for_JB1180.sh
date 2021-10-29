# This script is used to assemble JB1180 with smrtpipe.

# Enter into smrtshell environment with username "smrtanalysis"
cd /home/suofang/Wild_Strain/Separate/1116/0Rawdata   # because this folder can be accessed by smrtanalysis
[suofang@localhost 0Rawdata]$ bash
[suofang@localhost 0Rawdata]$ SMRT_ROOT=/opt/smrtanalysis
[suofang@localhost 0Rawdata]$ source $SMRT_ROOT/current/etc/setup.sh
(smrtanalysis-2.3.0) [suofang@localhost 0Rawdata]$ $SMRT_ROOT/smrtcmds/bin/smrtshell
(smrtshell-2.3.0) [smrtanalysis@localhost ~]$ SMRT_ROOT=/opt/smrtanalysis/




# construct JB1180 fofn file "JB1180_pacbio.fofn" shown below
# /home/suofang/Wild_Strain/Separate/1116/0Rawdata/m170807_070522_42278_c101236212550000001823297112091747_s1_p0.1.bax.h5
# /home/suofang/Wild_Strain/Separate/1116/0Rawdata/m170807_070522_42278_c101236212550000001823297112091747_s1_p0.2.bax.h5
# /home/suofang/Wild_Strain/Separate/1116/0Rawdata/m170807_070522_42278_c101236212550000001823297112091747_s1_p0.3.bax.h5

# transform input file
(smrtshell-2.3.0) [suofang@localhost 0Rawdata]$ fofnToSmrtpipeInput.py JB1180_pacbio.fofn > JB1180_input.xml

# do assembly, all parameters are stored in JB1180_hgap3_params.xml (find it in the filefolder)
# change options:  minLongReadLength=1000; genomeSize=12500000
(smrtshell-2.3.0) [suofang@localhost 0Rawdata]$ smrtpipe.py --params=JB1180_hgap3_params.xml xml:JB1180_input.xml


