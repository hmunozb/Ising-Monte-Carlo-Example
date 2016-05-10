Some Monte Carlo simulation code for the 2D ising model, tactically compiled on different machines with different jobs on the main function so that it doesn't take forever to simulate.

Output from this code was analyzed for my final report in PHYS 327 Experimental Physics at TAMU.

This code utilizes:
    - Concurrent Queue, Copyright (c) 2013-2016, Cameron Desrochers. Simplified BSD license in concurrentqueue.h
        http://moodycamel.com/blog/2014/a-fast-general-purpose-lock-free-queue-for-c++
    - Random-Number Utilities (randutil), Copyright (c) 2015 Melissa E. O'Neill.  The MIT license in randutils.hpp
        http://www.pcg-random.org/posts/developing-a-seed_seq-alternative.html
    - Easylogging++,  Copyright (c) 2015 muflihun.com 
        https://github.com/easylogging/easyloggingpp
    - The Lean Mean C++ Option Parser, Copyright (C) 2012 Matthias S. Benkmann
