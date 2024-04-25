#ifndef BETA_KMEANS_H
#define BETA_KMEANS_H

/* Authors: Greg Hamerly and Jonathan Drake
 * Feedback: hamerly@cs.baylor.edu
 * See: http://cs.baylor.edu/~hamerly/software/kmeans.php
 * Copyright 2014
 *
 * NaiveKmeans is the standard k-means algorithm that has no acceleration
 * applied. Also known as Lloyd's algorithm.
 */

#include "original_space_kmeans.h"

class BetaKmeans : public OriginalSpaceKmeans {
    public:
        virtual std::string getName() const { return "beta"; }
        virtual ~BetaKmeans() { free(); }
    protected:
        virtual int runThread(int threadId, int maxIterations);
        double *pointNorms;
        int *bins;
        int *beta;
        int *alpha;
        double *probabilities;
};

#endif

