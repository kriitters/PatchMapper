# PatchMapper
A C tool to identify, map, and analyze contiguous patches on 8-bit raster images.

## Contents
PatchMapper.c  compile on linux only, see Guide for gcc command  
PatchMapperGuide.pdf  
Example run:  
- inputmap - a 4177x6448 raster image, 8bit unsigned int, value 255 is "missing" (required)  
- parfile.txt, PatchMapper parameter file (required)  
- size.txt, specifies nrows and ncols in inputmap (required)  
- recode.txt, specifies how inputmap values should be re-coded (optional)  
