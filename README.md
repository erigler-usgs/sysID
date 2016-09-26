# sysID
E. J. Rigler's old Octave programs for system identification

## Summary

This is a collection of Octave (https://www.gnu.org/software/octave/) code I started
writing in graduate school, and used for my thesis and various radiation belt-related
research projects until ~2006. All the code here is intended for time-domain time 
series analysis and forecasting. It includes a variety of multi-input, multi-output, 
adaptive linear filters. However, I took several detours in my space weather research 
career, first to the ionosphere, then the sun, and finally down to the Earth's surface 
where I currently study geomagnetism at the USGS' Geological Hazards Science Center.

I've wanted to resurrect this code base for some time now, clean it up, apply minor 
changes necessary to make it Matlab(TM) compatible, and possibly translate it to Python,
my preferred numerical analysis environment these days. For now at least, all I can
manage is to convert my old RCS repository into a more modern git repo, and post it
Online at GitHub.

In addition to needing a thorough software review and overhaul (hey, some of this stuff
is over 15 years old!), this repository is in desparate need of documentation. Until
then, I've included my thesis. The "science" in the thesis is probably not all that 
interesting, but you might find useful information and context for the code herein
in Chapter 3, as well as Appendix C. In addition, you might look at the following
papers I wrote back near the beginning of my research career:

http://onlinelibrary.wiley.com/doi/10.1029/2003SW000036/abstract

http://onlinelibrary.wiley.com/doi/10.1029/2006JA012181/abstract

http://www.sciencedirect.com/science/article/pii/S1364682608000369

If you notice something really wrong with any of it, please describe in detail in a 
pull-request, and I will try to apply the fix. Until such time as I can make this 
project a priority again however, this repository should be considered static.

Happy Coding,
Josh Rigler
