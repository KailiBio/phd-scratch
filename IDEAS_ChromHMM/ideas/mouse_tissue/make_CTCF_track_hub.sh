#!/bin/bash

# -- Kaili
# This script is for getting CTCF peaks, then make track hub.

scriptDir="/data/zusers/fankaili/github/weng-lab/Kaili/IDEAS_ChromHMM/ideas/mouse_tissue/"
workDir="/data/zusers/fankaili/ideas/"

cd ${workDir}

# 1. get bigWig file
grep "CTCF" mm10_tissue_used_list.txt | grep "_14.5_day" > mm10_tissue_used_list_CTCF_peak_list.txt
grep "CTCF" mm10_tissue_used_list.txt | grep "_0_day" >> mm10_tissue_used_list_CTCF_peak_list.txt
#
vim mm10_tissue_used_list_CTCF_peak_list.txt
# manually add bigBed fileID

# 2. make track file

# hub file
if [ -f /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt ]
then
    rm /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
    touch /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
else
    touch /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
fi
#
echo "hub e14.5p0_CTCF_peak" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
echo "shortLabel e14.5p0_CTCF_peak" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
echo "longLabel e14.5p0_CTCF_peak" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
echo "genomesFile genomes_CTCF_peak_hub.txt" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
echo "email fankaili.bio@gmail.com" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt


# genome file
if [ -f /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/genomes_e14.5p0_CTCF_peak.txt ]
then
    rm /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/genomes_e14.5p0_CTCF_peak.txt
    touch /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/genomes_e14.5p0_CTCF_peak.txt
else
    touch /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/genomes_e14.5p0_CTCF_peak.txt
fi
#
echo "genome mm10" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/genomes_e14.5p0_CTCF_peak.txt
echo "trackDB trackDb_CTCF_peak_hub.txt" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/genomes_e14.5p0_CTCF_peak.txt

# trackDb file
mkdir CTCF_peak_hub
cd /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub
#
if [ -f /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt ]
then
    rm /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    touch /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
else
    touch /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
fi
#
while read line
do
    biosample=`awk '{print $1}' <<< ${line}`
    fileID=`awk '{print $6}' <<< ${line}`
    echo ${fileID}
    # download file
    wget https://www.encodeproject.org/files/${fileID}/@@download/${fileID}.bigBed
    # make trackDb file
    echo "track "${biosample} >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "priority 100" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "type bigBed" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    #echo "itemRgb on" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "maxItems 100000" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "bigDataUrl https://users.wenglab.org/fankaili/IDEAS/CTCF_peak_hub/${fileID}.bigBed" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "shortLabel "${biosample} >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "longLabel "${biosample}"_CTCF_peak" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "color 200,0,250" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "visibility dense" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
    echo "" >> /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/trackDb_e14.5p0_CTCF_peak.txt
done < ${workDir}mm10_tissue_used_list_CTCF_peak_list.txt



# 3. connect
# https://users.wenglab.org/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt
