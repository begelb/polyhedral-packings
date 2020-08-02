# polyhedral-packings
This is the code used to generate data on polyhedral packings that was created during the 2020 DIMACS REU. 

How to use this code:

Generate 3-connected planar graph data from plantri here: http://combos.org/plantri 

Save the output as a textfile named "PLANTRI_FORMAT.txt"

IN THE SAGE FILE: 

If you change the name of "PLANTRI_FORMAT.TXT', you must change the filename on line 3 to match.

Change "MATHEMATICA_FORMAT.txt" to whatever you want to name the file outputted by this code. You can then use the mathematica program ("Template_Math") to read that file.

This code is set up so that is convenient to handle many graphs at a time. In particular, Template_Math is set up to handle 250 graphs at a time. This is done using the function "Do[***lots of code***, 250]. Change the number of graphs by changing the number 250 on the last line. In the Sage code, Line 1 in the third box controls which graphs from PLANTRI_FORMAT.txt are being handled. For example, if you generate 1000 circle packings from 10-vertex graphs, we recommend doing it this way:

Use http://combos.org/plantri to output a text file with the plantri data (rename it PLANTRI_FORMAT.txt).

Change Line 1 of the third box to "for i in range(1, 251):"

Change "MATHEMATICA_FORMAT.txt" to "10_250.txt"

Use replace all in Template_Math in Mathematica so that Mathematica is reading the correct file in the ReadLine function used on the first 20 or so lines. You will need to specify the full path to 10_250.txt on your computer. Google "copy full path"+your type of computer to find out how to do that.

Run the modified Template_Math

Then repeat by changing Line 1 of the third box to "for i in range(250, 501):"

Change "MATHEMATICA_FORMAT.txt" to "10_500.txt" ...

If you have any questions I can be reached at begelb at muhlenberg dot edu
