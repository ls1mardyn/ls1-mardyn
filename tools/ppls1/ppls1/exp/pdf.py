# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 12:16:12 2019

@author: mheinen
"""

import re
import PyPDF2
import os
import numpy as np

def exp_pdf_merge(inppath,pattern,outpath,outfname):
    pdf_writer=PyPDF2.PdfFileWriter()
    flist=[fname for fname in os.listdir(inppath) if bool(re.search(pattern,fname))]
    flist=sorted(np.array(flist))
	
    for file in flist:
        fpath=os.path.join(inppath,file)
        pdf_reader=PyPDF2.PdfFileReader(fpath)
        #Opening each page of the PDF
        for pagenum in range(pdf_reader.numPages):
            page=pdf_reader.getPage(pagenum)
            pdf_writer.addPage(page)
    with open(os.path.join(outpath,outfname),'wb') as out:
        pdf_writer.write(out)
    print("Saved merged pdf as: "+os.path.join(outpath,outfname))
    return flist
