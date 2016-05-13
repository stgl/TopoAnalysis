#!/usr/bin/env python

Ao = 250000
theta = 0.4

import demMethods as d

d.processForTheta('af', Ao, theta)
d.processForTheta('as', Ao, theta)
d.processForTheta('au', Ao, theta)
d.processForTheta('ca', Ao, theta)
d.processForTheta('eu', Ao, theta)
d.processForTheta('na', Ao, theta)
d.processForTheta('sa', Ao, theta)

theta = 0.6

d.processForTheta('af', Ao, theta)
d.processForTheta('as', Ao, theta)
d.processForTheta('au', Ao, theta)
d.processForTheta('ca', Ao, theta)
d.processForTheta('eu', Ao, theta)
d.processForTheta('na', Ao, theta)
d.processForTheta('sa', Ao, theta)

theta = 0.5

d.processForTheta('as', Ao, theta)
