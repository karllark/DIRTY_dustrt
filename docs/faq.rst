##########################
Frequently Asked Questions
##########################

General
=======

What does DIRTY stand for?
--------------------------

DIRTY stands for DustI Radiative Transfer, Yeah!.
The name was hatched at Chimes, a pub in Baton Rouge, Louisiana.

What is DIRTY?
--------------

DIRTY is a Monte Carlo dust radiative transfer/emission code that computes the
radiative transfer and dust emission from arbitrary distributions of dust
illuminated by arbitrary distributions of sources (usually stars).

Who wrote/writes DIRTY?
-----------------------

Mainly Karl Gordon (radiative transfer, optical dust emission, I/O) and
Karl Misselt (dust grain model, thermal/IR dust emission).
Additional coding help (fixing the Karl's errors) has been given by Ka-Hei Law.
Many others have helped by using the model and giving feedback.

Why C++?
--------

DIRTY is a computationally intensive code and is written in C++.  Speed is
very much needed.  Most (if not all) similar models are written in a
language that can be compiled.  Scripting languages like python and IDL
are not well suited to such models.

DIRTY support?
--------------

No assumption of support should be made.  Support is provide on the whim
of the Karl's.  The Karl's are interested in research using DIRTY
and interesting (to them) research can get support.

Where are the algorithms used by DIRTY described?
-------------------------------------------------

General organization and code:
`Gordon et al. (2001, ApJ, 551, 269)
<https://ui.adsabs.harvard.edu/#abs/2001ApJ...551..269G/abstract>`_,
`Misselt et al. (2001, ApJ, 551, 277)
<https://ui.adsabs.harvard.edu/#abs/2001ApJ...551..277M/abstract>`_

Some updates:
`Law et al. (2018, ApJS, 236, 32)
<https://ui.adsabs.harvard.edu/#abs/2018ApJS..236...32L/abstract>`_

Getting Started
===============

Where do I start if I want to run the DIRTY model?
--------------------------------------------------

See :ref:`install`.

Where do I find documentation on setting up the DIRTY parameter file?
---------------------------------------------------------------------

Check out the (not yet complete) documentation at DirtyParameterFile.
Feel free to add to the documentation if it is unclear. This would be greatly appreciated.

What do I do when DIRTY crashes/has an error?
---------------------------------------------

Open an issue on the `GitHub issue tracker
<https://github.com/karllark/DIRTY_dustrt/issues>`_ with the
parameter file and the error output from the DIRTY code (cut and paste from the terminal).
Also include the version of the code you are using and any other information
that you think might help us debug the problem.

Testing
=======

Has DIRTY been benchmarked?
---------------------------

The DIRTY code has been benchmarked over the years against various
1D and 2D cases.  More recently, DIRTY participated in the 1st ever
fully 3D dust radiative transfer benchmark and produced very similar
results to the other 6 codes
(`Gordon et al. 2017, A&A, 603A, 114
<https://ui.adsabs.harvard.edu/#abs/2017A&A...603A.114G/abstract>`_).
