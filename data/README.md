Downloading bam-readcount data from DNAnexus

ran commands from inside directories

_RM data_
```
cat DNAnexus_export_urls-20250402-162218.txt| grep "snv" | parallel -j 12 -d "\r\n" wget
```

_nonRM data_
```
cat DNAnexus_export_urls-20250402-162552.txt| grep "snv" | parallel -j12 -d "\r\n" wget
``` 
