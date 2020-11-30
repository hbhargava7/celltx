# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell A. Lim
# University of California, San Francisco

import time
import numbers
import math


def format_timedelta(diff):
    if isinstance(diff, numbers.Real):
        s = diff
        d = 0
    else:
        s = diff.seconds
        d = diff.days
    h,r = divmod(s, 3600)
    m,s = divmod(r, 60)
    if d > 0:
        return '%02d:%02d:%02d:%02ds' % (diff.days, h,m,s)
    else:
        return '%02d:%02d:%02ds' % (h,m,s)


