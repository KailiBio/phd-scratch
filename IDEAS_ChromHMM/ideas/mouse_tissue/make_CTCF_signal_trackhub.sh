#!/bin/bash

# -- Kaili
# This script is for making CTCF signal track hub.

mkdir /data/public_html_users/fankaili/IDEAS/CTCF_signal_hub/
cd /data/public_html_users/fankaili/IDEAS/CTCF_signal_hub/

# .hub file
cp /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/hub_e14.5p0_CTCF_peak.txt hub_e14.5p0_CTCF_signal.txt
vim hub_e14.5p0_CTCF_signal.txt

# .genome file
cp /data/public_html_users/fankaili/IDEAS/CTCF_peak_hub/genomes_e14.5p0_CTCF_peak.txt genomes_e14.5p0_CTCF_signal.txt
vim genomes_e14.5p0_CTCF_signal.txt

#
if [ -f trackDb_e14.5p0_CTCF_signal.txt ]
then
    rm trackDb_e14.5p0_CTCF_signal.txt
    touch trackDb_e14.5p0_CTCF_signal.txt
else
    touch trackDb_e14.5p0_CTCF_signal.txt
fi
#
while read line
do
    sample=`awk '{print $1}' <<< ${line}`
    expID=`awk '{print $4}' <<< ${line}`
    fileID=`awk '{print $5}' <<< ${line}`
    # get .bigwig file
    cp /data/projects/encode/data/${expID}/${fileID}.bigWig ./
    #
    echo "track "${sample} >> trackDb_e14.5p0_CTCF_signal.txt
    echo "shortLabel "${sample} >> trackDb_e14.5p0_CTCF_signal.txt
    echo "longLabel "${sample} >> trackDb_e14.5p0_CTCF_signal.txt
    echo "priority 200" >> trackDb_e14.5p0_CTCF_signal.txt
    echo "type bigWig" >> trackDb_e14.5p0_CTCF_signal.txt
    echo "maxItems 100000" >> trackDb_e14.5p0_CTCF_signal.txt
    echo "bigDataUrl https://users.wenglab.org/fankaili/IDEAS/CTCF_signal_hub/"${fileID}".bigWig" >> trackDb_e14.5p0_CTCF_signal.txt
    echo "color 200,0,250" >> trackDb_e14.5p0_CTCF_signal.txt
    echo "visibility full" >> trackDb_e14.5p0_CTCF_signal.txt
    echo "autoScale on" >> trackDb_e14.5p0_CTCF_signal.txt
    echo "" >> trackDb_e14.5p0_CTCF_signal.txt
done < /data/zusers/fankaili/ideas/mm10_tissue_used_list_CTCF_peak_list.txt
