import sys, os
import attracthtmlform 
import spyder
from attractmodel import AttractEasyModel
import formeasy as form

cgi = sys.argv[1]

f = AttractEasyModel._form()
f = form.webform(f)
html = attracthtmlform.htmlform(
 form=f, cgi=cgi, 
 header=form.header, footer=form.footer, header_indentation = 12, 
 newtab=True
)
print(html)