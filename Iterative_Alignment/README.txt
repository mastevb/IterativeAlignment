This doc will help you run the iterative alignment program. This program essentially changes a similarity matrix to recognize when strong binding peps are more similar to eachother and weaker binders are less similar to strong binders.

For each dataset, you must first obtain the top 10 and bottom 10 percent of the binders from each list, and label these high.csv and low.csv. Place these files in a file labeled "1" under the peps folder. 


Then run the file from the command console using "python IterativeAlignment.py" after navigating to the folder.

The windows will ask you to locate the saves, peps and matrix folders. 

Then it will ask you for some numbers. 

For iterations, use 10000. (How many rounds of modification will occur?)
For cutoff, use 5000. (How many incorrect changes can be made before the matrix is considered trained?)
For dataset selection, use 10000. (This asks how much of the data from the high/low lists to use, this is an older unneccessary feature.)
For perturb decimal, use .5 (How much should a particular value be changed during each round?)

When you run it, two numbers should be increasing in the command window. One should be resetting every now and then. If it reaches 5000, the program will stop and the matrix will be considered trained.