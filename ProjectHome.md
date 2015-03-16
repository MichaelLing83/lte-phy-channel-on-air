# Important Note #
This project has been split into two:
  1. The visualizer part is moved to [pyLTEVisualizer](http://code.google.com/p/pyltevisualizer/)
  1. The simulation according to 3GPP part is moved to [pyLTESimulator](http://code.google.com/p/pyltesimulator/)


---


Why did I create this project?

Well, at the beginning of learning TD-LTE technology, I really had problem visualizing how the air interface of TDD-LTE system looks like.  So I tried to depict these channels. And this project is the work I've done.

What is the status of this project?

In July, 2010, I left TD-LTE area and started working in WCDMA, so this project was suspended from then. This year, I'm back in TD-LTE again, so it gets revived! I'll keep updating this project, as long as I'm in this area, and thank you for your time to check it out.

How can I use this project?

Just download the featured package on project home page, and check out the image files. If you want to change configurations and run the Python codes, you have to install Python and  the modules mentioned below. I would like to recommend you to install Python(x,y) from http://code.google.com/p/pythonxy/ and LEO from http://webpages.charter.net/edreamleo/front.html. Python(x,y) is a little large but it has the whole package.

What is the goal of this project?

My first goal is to write a visualizer to depict various DL/UL signals and channels in the time-frequency resource grid. And after getting deeper in this technology, I now have a more ambitious goal that is to simulate various signals and channels. See how they are sent and how they can be detected and used. This could be too ambitious, but I want to see how far I can go on this.

How is the source code organized?

This project is based on Python, uses PIL to output resource grid image, uses Numpy and SciPy to do signal processing and scientific calculation, and uses Matplotlib to draw diagrams of simulation results. The code is organized in LEO, a great literature programming tool written in Python. So if you can install and use LEO to open the lte-phy-channel-on-air.leo file, and all source codes are there. If you don't want to use LEO, I think you can just read `*`.py files directly, but without LEO the source codes are quite messy.