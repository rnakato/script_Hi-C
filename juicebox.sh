#!/bin/bash

singularity exec rnakato_juicer.img java -Xms512m -Xmx2048m -jar /opt/Juicebox.jar
#singularity exec rnakato_juicer.img java -Xms512m -Xmx2048m -jar /opt/juicer_tools_1.14.08.jar hiccups -r 5000,10000,25000 -k KR  /home/rnakato/Hi-C/Sakata_RPE/JuicerResults/CTCFKD_1/aligned/inter_30.hic test
