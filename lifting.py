#Lifting Progression

from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

# def date_nums(date_list):
#     date_nums = []
#     current_year = 2024
#     for date_str in date_list:
#         date_object = datetime.strptime(date_str, "%b %d %Y")
#         if date_object.year==2022:
#             day_of_year = date_object.timetuple().tm_yday
#         if date_object.year==2023:
#             day_of_year = 365+date_object.timetuple().tm_yday
#         if date_object.year==2024:
#             day_of_year = 365*2+date_object.timetuple().tm_yday
#         date_nums.append(day_of_year)
#     return date_nums


def date_nums(date_list):
    date_nums = []
    current_year = 2024
    for date_str in date_list:
        date_object = datetime.strptime(date_str, "%b %d %Y")
        if date_object.year==2015:
            day_of_year = 365*0+date_object.timetuple().tm_yday
        if date_object.year==2016:
            day_of_year = 365*1+date_object.timetuple().tm_yday
        if date_object.year==2017:
            day_of_year = 365*2+date_object.timetuple().tm_yday
        if date_object.year==2018:
            day_of_year = 365*3+date_object.timetuple().tm_yday
        if date_object.year==2019:
            day_of_year = 365*4+date_object.timetuple().tm_yday
        if date_object.year==2020:
            day_of_year = 365*5+date_object.timetuple().tm_yday
        if date_object.year==2021:
            day_of_year = 365*6+date_object.timetuple().tm_yday
        if date_object.year==2022:
            day_of_year = 365*7+date_object.timetuple().tm_yday
        if date_object.year==2023:
            day_of_year = 365*8+date_object.timetuple().tm_yday
        if date_object.year==2024:
            day_of_year = 365*9+date_object.timetuple().tm_yday
        if date_object.year==2025:
            day_of_year = 365*10+date_object.timetuple().tm_yday
        if date_object.year==2026:
            day_of_year = 365*11+date_object.timetuple().tm_yday
        date_nums.append(day_of_year)
    return date_nums

# f, ax = plt.subplots(1, 2, figsize=(6, 4))
f, ax = plt.subplots(5, 1, figsize=(12, 9))
# f, ax = plt.subplots(1, 4, figsize=(16, 9))
plt.subplots_adjust(hspace=0.5)

squats={
    "Jul 21 2015":5*5*170,
    "Jul 25 2015":5*5*175,
    "Jul 27 2015":5*5*175,
    "Jul 28 2015":5*5*175,
    "Jul 30 2015":5*5*185,
    "Jul 31 2015":5*5*200, #probably wasn't 5x5
    "Aug 1 2015":5*5*190,
    "Aug 11 2015":5*5*175,
    "Aug 12 2015":5*5*180,
    "Aug 22 2015":5*5*180,
    "Sep 1 2015":5*5*180,
    "Sep 2 2015":5*5*185,
    "Nov 23 2017":4*5*185,
    "Oct 19 2022":(4*4*185),
    "Oct 21 2022":(5*3*185+3*185),
    "Nov 16 2022":(5*135+5*185+5*195+5*205+4*205),
    "Dec 2 2022":(10*2*135),
    "Dec 4 2022":(10*2*135),
    "Dec 6 2022":(5*5*205),
    "Dec 8 2022":(5*5*210),
    "Dec 12 2022":(5*215+5*220+5*225+4*225),
    "Dec 16 2022":(5*2*225+4*225),
    "Dec 18 2022":(5*3*225+3*225),
    "Dec 20 2022":(5*2*225),
    "Jan 18 2023":(5*2*185+5*1*205),
    "Apr 30 2023":(5*1*185+2*1*225+1*225+5*2*185),
    "Jul 24 2023":(5*5*185),
    "Jul 27 2023":(5*4*205),
    "Dec 30 2023":(5*2*205+4*205+3*205),
    "Jan 23 2024":(5*2*205),
    "Jan 24 2024":(5*3*185),
    "Jan 25 2024":(5*4*205),
    "Feb 2 2024":(5*5*205),
    "Feb 9 2024":(5*1*205+5*1*215+5*1*220+2*1*220),
    "Feb 15 2024":(5*205+5*2*215+5*205),
    "Feb 21 2024":(5*2*215+5*3*225),
    "Feb 26 2024":(5*5*225),
    "Feb 29 2024":5*5*230,
    "Mar 5 2024":5*5*230,
    "Mar 8 2024":5*3*235,
    "Mar 20 2024":5*3*235+5*2*225,
    "Apr 5 2024":5*5*225,
    "Apr 10 2024":5*5*225,
    "Apr 16 2024":5*5*235,
    "Apr 22 2024":5*5*240,
    "Aug 22 2024":5*1*205+5*3*225+3*1*225,
    "Aug 29 2024":5*5*225,
    "Sep 18 2024":5*5*225,
    "Oct 7 2024":5*2*225+5*3*135,
    "Mar 10 2025":5*5*205,
    "Nov 11 2025":5*1*135+5*2*185,
    "Dec 8 2025":5*4*205+4*1*205,
    "Dec 15 2025":5*5*215,
    "Dec 21 2025":5*5*225,
    "Jan 6 2026":5*5*230,
    "Jan 12 2026":5*5*235,
    "Feb 2 2026":5*5*240,
    "Feb 10 2026":5*5*235,
    "Feb 16 2026":5*5*240,
    }
squats5x5={
    "Feb 2 2024":(5*5*205)/25,
    "Feb 26 2024":(5*5*225)/25,
    "Feb 29 2024":5*5*230/25,
    "Mar 5 2024":5*5*230/25,
    "Apr 5 2024":5*5*225/25,
    "Apr 10 2024":5*5*225/25,
    "Apr 16 2024":5*5*235/25,
    "Apr 22 2024":5*5*240/25,
    "Aug 29 2024":5*5*225/25,
    "Sep 18 2024":5*5*225/25,
    "Mar 10 2025":5*5*205/25,
    "Dec 08 2025":5*5*205/25, #actually only did 4 reps on last set, but including for granularity
    "Dec 15 2025":5*5*215/25,
    "Dec 21 2025":5*5*225/25,
    "Jan 6 2026":5*5*230/25,
    "Jan 12 2026":5*5*235/25,
    "Feb 2 2026":5*5*240/25,
    "Feb 10 2026":5*5*235/25,
    "Feb 16 2026":5*5*240/25,
    }

bench={
    "Jul 25 2015":5*5*135, #probably wasn't 5x5
    "Jul 28 2015":5*5*140, #probably wasn't 5x5
    "Jul 31 2015":5*5*140, #probably wasn't 5x5
    "Aug 11 2015":5*5*120, 
    "Aug 22 2015":5*5*120,
    "Sep 1 2015":5*5*135,
    "Oct 19 2022":(4*4*135),
    "Dec 2 2022":(10*2*115),
    "Dec 6 2022":(5*135+5*165+5*155+2*155+3*135),
    "Dec 12 2022":(4*165+3*175+2*175+4*155+3*155),
    "Dec 18 2022":(3*165+5*2*165+3*165+3*135),
    "Apr 30 2023":(5*135+4*165+4*155+4*135),
    "Jul 24 2023":(5*4*135),
    "Dec 30 2023":(5*4*135+3*135),
    "Jan 23 2024":(5*1*135+5*2*155),
    "Jan 25 2024":(5*1*135+4*2*155+5*1*140),
    "Jan 30 2024":(5*1*155+4*1*155+5*1*155+2*1*155+4*1*135),
    "Feb 7 2024":(5*5*135),
    "Feb 13 2024":(5*5*155),
    "Feb 16 2024":(5*2*160+4*160+3*2*155),
    "Feb 22 2024":(5*4*160+4*160),
    "Feb 27 2024":(5*2*165+4*165+4*155+3*155),
    "Mar 1 2024":5*155+5*160+4*3*160,
    "Mar 6 2024":5*5*165,
    "Mar 17 2024":5*155+5*165+5*3*170,
    "Mar 27 2024":5*2*170+4*2*170+3*170,
    "Apr 8 2024":5*155+5*3*160+4*160,
    "Apr 11 2024":5*2*170+5*165+4*165+3*165,
    "Apr 17 2024":5*5*165,
    "Apr 26 2024":5*5*170,
    "May 29 2024":5*1*165+5*2*160+4*2*155,
    "Jun 20 2024":5*1*135+5*4*155,
    "Aug 19 2024":5*4*135,
    "Aug 23 2024":5*4*155,
    "Sep 3 2024":5*5*160,
    "Sep 10 2024":5*1*135+5*4*155,
    "Sep 19 2024":5*5*165,
    "Sep 30 2024":5*5*165,
    "Nov 8 2024":5*5*135,
    "Feb 24 2025":5*1*135+5*1*145+5*1*150+5*1*155+4*1*155,
    "Feb 28 2025":5*1*135+2*1*155+5*2*135+3*1*135,
    "Mar 6 2025":5*1*135+5*2*155+4*2*155,
    "Mar 11 2025":5*1*135+5*1*145+4*1*145+3*2*145,
    "Mar 21 2025":5*5*145,
    "Mar 24 2025":5*1*135+5*4*150,
    "Mar 31 2025":5*4*155+3*1*155,
    "Apr 8 2025":5*1*135+5+1+155+5*3*145,
    "Apr 12 2025":5*5*155,
    "Apr 23 2025":5*5*160,
    "Apr 26 2025":5*2*160+5*2*155+4*1*155,
    "Apr 30 2025":5*5*160,
    "May 4 2025":5*4*165+4*1*165,
    "May 7 2025":5*5*165,
    "May 14 2025":5*5*165,
    "May 28 2025":5*2*165+5*3*155,
    "Jun 4 2025":5*5*165,
    "Jun 12 2025":5*1*135+5*4*155,
    "Jun 18 2025":5*2*165+5*1*170+5*2*175,
    "Jun 24 2025":5*5*165,
    "Jun 30 2025":5*5*170,
    "Jul 8 2025":5*1*170+5*4*180,
    "Jul 12 2025":5*1*135+5*1*180+5*1*175+5*1*155,
    "Jul 23 2025":5*5*175,
    "Jul 29 2025":5*3*175+3*2*175,
    "Aug 5 2025":5*5*170,
    "Aug 19 2025":2*1*175+5*2*155+4*2*155,
    "Aug 27 2025":5*4*170+4*1*170,
    "Aug 29 2025":5*3*175+3*1*175+5*1*155,
    "Sep 1 2025":5*1*135+5*2*165+4*1*165,
    "Sep 16 2025":5*3*175+3*1*165+5*1*155,
    "Nov 6 2025":5*3*165+4*1*165+3*1*165,
    "Nov 21 2025":5*1*135+5*1*155+5*3*165,
    "Dec 4 2025":5*2*155+5*3*160,
    "Dec 12 2025":5*4*160+3*1*160,
    "Dec 20 2025":5*5*165,
    "Dec 29 2025":5*5*175,
    "Jan 9 2026":5*5*180,
    "Jan 27 2026":5*5*185,
    "Jan 30 2026":5*5*190,
    "Feb 06 2026":5*5*185,
    "Feb 14 2026":5*5*195,
    }
bench5x5={
    "Feb 7 2024":(5*5*135)/25,
    "Feb 13 2024":(5*5*155)/25,
    "Mar 6 2024":5*5*165/25,
    "Apr 17 2024":5*5*165/25,
    "Apr 26 2024":5*5*170/25,
    "Sep 3 2024":5*5*160/25,
    "Sep 19 2024":5*5*165/25,
    "Sep 30 2024":5*5*165/25,
    "Nov 8 2024":5*5*135/25,
    "Mar 21 2025":5*5*145/25,
    "Apr 12 2025":5*5*155/25,
    "Apr 23 2025":5*5*160/25,
    "Apr 30 2025":5*5*160/25,
    "May 7 2025":5*5*165/25,
    "May 14 2025":5*5*165/25,
    "Jun 4 2025":5*5*165/25,
    "Jun 24 2025":5*5*165/25,
    "Jun 30 2025":5*5*170/25,
    "Jul 23 2025":5*5*175/25,
    "Aug 5 2025":5*5*170/25,
    "Dec 4 2025":5*5*155/25, #slightly fake number, but included for granularity
    "Dec 20 2025":5*5*165/25,
    "Dec 29 2025":5*5*175/25,
    "Jan 9 2026":5*5*180/25,
    "Jan 27 2026":5*5*185/25,
    "Jan 30 2026":5*5*190/25,
    "Feb 06 2026":5*5*185/25,
    "Feb 14 2026":5*5*195/25,
    "Feb 21 2026":5*5*195/25,
    }

ohp={
    "Jul 25 2015":5*5*75,
    "Jul 28 2015":5*5*75,
    "Jul 31 2015":5*5*80,
    "Aug 11 2015":5*5*75,
    "Aug 22 2015":5*5*75,
    "Sep 1 2015":5*5*75,
    "Nov 23 2017":4*5*95,
    "Oct 21 2022":(3*95+4*95+2*95),
    "Dec 2 2022":(10*2*75),
    "Dec 4 2022":(10*2*75),
    "Dec 8 2022":(5*3*95+4*95),
    "Dec 16 2022":(3*95+5*95+4*95+3*2*95),
    "Dec 20 2022":(3*95+5*95+4*2*95),
    "Jan 18 2023":(5*1*95+4*2*95),
    "Jul 27 2023":(3*95+4*95+3*4*95),
    "Dec 30 2023":(95+3*2*85+4*85),
    "Jan 23 2024":(6*2*75+3*1*75),
    "Jan 26 2024":(3*4*95),
    "Feb 13 2024":(3*5*85),
    "Feb 16 2024":(5*85+3*3*85+5*65),
    "Feb 22 2024":(4*2*85+3*80+4*2*75),
    "Feb 27 2024":(4*80+5*4*75),
    "Mar 1 2024":5*85+5*90+5*3*95,
    "Mar 4 2024": 5*85+5*90+4*90+5*90+4*90,
    "Mar 17 2024":5*85+5*4*95,
    "Mar 27 2024":3*95+95+5*2*85+3*2*85,
    "Apr 8 2024":4*5*95, #did OHP before bench
    "Apr 11 2024":5*2*95+3*95+5*65+5*65, #did OHP after bench
    "Apr 17 2024":5*5*95, #did OHP before bench
    "Apr 26 2024":4*2*95+3*3*95, #did OHP after bench
    "May 29 2024":5*4*95+3*1*95, #did OHP before bench
    "Jun 20 2024":5*5*85,
    "Aug 19 2024":5*2*75,
    "Aug 23 2024":5*3*85+4*1*85+5*1*85,
    "Sep 3 2024":5*3*85+4*2*85,
    "Sep 10 2024":5*5*85,
    "Sep 19 2024":5*5*85,
    "Sep 30 2024":5*3*90+4*2*90,
    "Nov 8 2024":5*3*85+4*2*85,
    "Feb 24 2025":5*3*85+4*1*85+3*1*85,
    "Feb 28 2025":5*2*85+4*1*85+3*2*85,
    "Mar 6 2025":3*1*85+5*2*75+4*1*75+3*1*75,
    "Mar 11 2025":5*5*75,
    "Mar 21 2025":5*3*80+4*1*80+3*1*80,
    "Mar 24 2025":5*4*80+4*1*80,
    "Mar 31 2025":5*4*80,
    "Apr 4 2025":5*5*80,
    "Apr 8 2025":5*5*75,
    "Apr 12 2025":5*4*85+4*1*85,
    "Apr 23 2025":5*3*85+4*2*85,
    "Apr 26 2025":5*5*90,
    "Apr 30 2025":5*4*95+4*1*95,
    "May 4 2025":5*5*95,
    "May 7 2025":5*1*95+5*3*100+4*1*100,
    "May 14 2025":4*1*100+5*4*95,
    "May 27 2025":4*1*105+5*4*95,
    "May 28 2025":5*2*95+5*2*85,
    "Jun 4 2025":5*1*105+5*1*95+3*1*95+5*1*95+4*1*95,
    "Jun 12 2025":5*4*95+4*1*95,
    "Jun 18 2025":5*3*100+4*2*95,
    "Jun 24 2025":5*5*95,
    "Jun 30 2025":5*2*100+4*3*100,
    "Jul 8 2025":4*1*100+5*3*95+4*1*95,
    "Jul 12 2025":5*4*105+3*1*105,
    "Jul 23 2025":5*2*95+4*3*95,
    "Jul 29 2025":3*1*105+5*4*95,
    "Aug 5 2025":5*1*105+4*1*105+3*2*105+1*1*105,
    "Aug 19 2025":5*2*105+3*3*105,
    "Aug 27 2025":4*3*95+3*2*95,
    "Sep 1 2025":5*3*100+4*2*100,
    "Sep 16 2025":3*1*90+5*1*90+3*2*90,
    "Oct 2 2025":4*4*95+3*1*95,
    "Nov 6 2025":5*4*85+3*1*85,
    "Nov 17 2025":5*3*85+5*2*95,
    "Nov 21 2025":3*4*95+5*1*95,
    "Dec 4 2025":3*1*95+5*3*85+4*1*85,
    "Dec 12 2025":5*4*95+4*1*95,
    "Dec 20 2025":5*3*95+3*1*95+2*1*95,
    "Dec 29 2025":4*2*95+3*2*95,
    "Jan 9 2026":5*5*90,
    "Jan 17 2026":5*5*95,
    "Jan 27 2026":5*5*95,
    "Jan 31 2026":5*5*100,
    "Feb 16 2026":5*5*105,
    }
ohp5x5={
    "Feb 27 2024":(5*5*75)/25, #doctored this number a bit
    "Mar 17 2024":5*5*85/25, #doctored this number a bit
    "Apr 17 2024":5*5*95/25, #did OHP before bench
    "Jun 20 2024":5*5*85/25,
    "Sep 10 2024":5*5*85/25,
    "Sep 19 2024":5*5*85/25,
    "Mar 11 2025":5*5*75/25,
    "Apr 4 2025":5*5*80/25,
    "Apr 8 2025":5*5*75/25,
    "Apr 26 2025":5*5*90/25,
    "May 4 2025":5*5*95/25,
    "Jun 24 2025":5*5*95/25,
    "Nov 6 2025":5*5*80/25, #slightly fake number but included for granularity 
    "Jan 9 2026":5*5*90/25,
    "Jan 17 2026":5*5*95/25,
    "Jan 27 2026":5*5*95/25,
    "Jan 31 2026":5*5*100/25,
    "Feb 16 2026":5*5*105/25,
    }

deadlifts={
    "Jul 21 2015":5*5*175,
    "Jul 27 2015":5*5*170,
    "Jul 30 2015":5*5*165,
    "Aug 1 2015":5*5*175,
    "Aug 12 2015":5*5*175,
    "Sep 2 2015":5*5*185,
    "Dec 6 2022":(3*3*205),
    "Dec 12 2022":(4*225+3*4*225),
    "Dec 18 2022":(5*235),
    "Apr 30 2023":(5*225+3*225),
    "Jul 24 2023":(5*225),
    "Jan 25 2024":(5*1*245),
    "Feb 2 2024":(5*1*255),
    "Feb 15 2024":(5*1*255),
    "Feb 23 2024":(5*225+5*4*245),
    "Feb 28 2024": 5*5*250,
    "Mar 4 2024":5*5*255,
    "Mar 7 2024":5*225+3*3*260+260,
    "Mar 28 2024":5*5*245,
    "Apr 9 2024":5*245+5*250+5*255+3*2*255,
    "Apr 12 2024":5*4*250+4*250,
    "Apr 19 2024":5*5*250,
    "May 30 2024":5*1*225+5*4*250,
    "Aug 20 2024":5*5*225,
    "Aug 28 2024":5*1*225+5*4*235,
    "Sep 11 2024":5*1*225+5*4*235,
    "Sep 24 2024":5*5*240,
    "Oct 2 2024":5*1*225+5*4*245,
    "Feb 25 2025":5*1*135+5*4*225,
    "Mar 4 2025":5*3*225,
    "Mar 20 2025":5*4*225+2*1*225,
    "Jun 22 2025":5*4*225,
    "Jun 27 2025":5*5*225,
    "Jul 2 2025":3*5*230, #shorter sets bc rear delt pain
    "Jul 10 2025":5*1*225+5*1*235+4*1*235+3*2*235, #used straight grip for every lift hereafter
    "Jul 17 2025":5*3*235, #started using 5x3 scheme instead of 5x5 for every lift hereafter
    "Aug 20 2025":5*3*225, #searing pain in left rear delt
    "Nov 8 2025":3*5*225,
    "Nov 19 2025":5*3*225,
    "Dec 3 2025":5*2*230+4*1*230,
    "Dec 10 2025":5*3*230,
    "Dec 17 2025":5*3*235,
    "Dec 28 2025":5*3*245,
    "Jan 8 2026":5*3*250,
    "Jan 14 2026":5*3*255,
    "Jan 22 2026": 5*3*265,
    "Jan 28 2026": 5*3*275,
    "Feb 4 2026": 5*3*280,
    }
deadlifts5x5={ #not exactly 5 sets for all and thats okay
    "Feb 28 2024": 5*5*250/25,
    "Mar 4 2024":5*5*255/25,
    "Mar 28 2024":5*5*245/25,
    "Apr 19 2024":5*5*250/25,
    "Aug 20 2024":5*5*225/25,
    "Sep 24 2024":5*5*240/25,
    "Mar 4 2025":5*3*225/15,
    "Jun 22 2025":5*4*225/20,
    "Jun 27 2025":5*5*225/25,
    "Jul 17 2025":5*3*235/15, #started using 5x3 scheme instead of 5x5 for every lift hereafter
    "Aug 20 2025":5*3*225/15, #searing pain in left rear delt
    "Nov 19 2025":5*3*225/15,
    "Dec 10 2025":5*3*230/15,
    "Dec 17 2025":5*3*235/15,
    "Dec 28 2025":5*3*245/15,
    "Jan 8 2026":5*3*250/15,
    "Jan 14 2026":5*3*255/15,
    "Jan 22 2026": 5*3*265/15,
    "Jan 28 2026": 5*3*275/15,
    "Feb 4 2026": 5*3*280/15,
    }

row={
    "Jul 21 2015":5*5*110,
    "Jul 27 2015":5*5*115,
    "Jul 30 2015":5*5*120,
    "Aug 1 2015":5*5*130,
    "Aug 12 2015":5*5*135,
    "Sep 2 2015":5*5*135,
    "Dec 4 2022":(10*2*115),
    "Dec 8 2022":(5*4*135),
    "Dec 16 2022":(5*5*135),
    "Dec 20 2022":(5*135+5*145+4*145+8*115+5*125),
    "Jul 27 2023":(5*5*135),
    "Jan 24 2024":(3*1*135+5*2*135+4*1*135),
    "Jan 26 2024":(5*2*115+5*2*120),
    "Jan 31 2024":(5*1*115+5*4*125),
    "Feb 8 2024":(5*4*135),
    "Feb 14 2024":(5*1*115+5*4*135),
    "Feb 20 2024":(5*3*135),
    "Feb 23 2024":(5*4*135+4*135),
    "Feb 28 2024": 5*2*135+5*145+4*145+5*135,
    "Mar 4 2024": 5*145+5*4*155,
    "Mar 7 2024":5*155+5*4*145,
    "Mar 19 2024":5*5*162,
    "Mar 28 2024":5*165+5*4*160,
    "Apr 9 2024":5*5*145,
    "Apr 19 2024":5*5*155,
    "May 30 2024":5*4*155,
    "Aug 20 2024":5*5*135,
    "Aug 28 2024":5*2*145+5*2*155+4*1*155,
    "Sep 11 2024":5*1*135+5*4*145,
    "Sep 24 2024":5*1*155+5*3*150,
    "Oct 2 2024":5*1*155+5*1*145+5*1*150+5*2*155,
    "Mar 4 2025":5*5*135,
    #"Mar 20 2025":8*5*95,
    "Mar 29 2025":5*4*145,
    "Apr 3 2025":5*5*145,
    "Apr 7 2025":5*5*150,
    "Apr 11 2025":5*1*135+5*4*155,
    "Apr 25 2025":5*5*155,
    "Apr 28 2025":5*4*160,
    "May 2 2025":5*5*160,
    "May 6 2025":5*5*165,
    "May 9 2025":5*1*170+5*4*165,
    "May 29 2025":5*5*165,
    "Jun 22 2025":5*1*135+5*2*155,
    "Jun 27 2025":5*1*155+5*4*145,
    "Jul 2 2025":5*5*155,
    "Jul 10 2025":3*1*165+5*2*160+5*1*155+3*1*155,
    "Jul 17 2025":5*5*155,
    "Jul 25 2025":5*5*155,
    "Aug 20 2025":5*5*135,
    "Aug 29 2025":5*2*135+5*1*145+5*2*150,
    "Sep 3 2025":5*3*145+5*2*150,
    "Sep 17 2025":8*5*135,
    "Nov 8 2025":10*3*115,
    "Nov 19 2025":5*5*135,
    "Dec 3 2025":5*5*140,
    "Dec 10 2025":5*5*145,
    "Dec 17 2025":5*5*150,
    "Dec 28 2025":5*5*155,
    "Jan 8 2026":5*4*160+4*1*160,
    "Jan 14 2026":5*5*160,
    "Jan 22 2026":5*5*165,
    "Jan 28 2026":5*5*170,
    "Feb 14 2026":5*5*170,
    "Feb 18 2026":5*5*175,
    }
row5x5={
    "Feb 14 2024":(5*5*130)/25, #doctored this number a bit
    "Mar 4 2024": 5*5*150/25, #doctored this number a bit
    "Mar 19 2024":5*5*162/25,
    "Mar 28 2024":5*5*160/25, #doctored this number a bit
    "Apr 9 2024":5*5*145/25,
    "Apr 19 2024":5*5*155/25,
    "Aug 20 2024":5*5*135/25,
    "Mar 4 2025":5*5*135/25,
    "Apr 3 2025":5*5*145/25,
    "Apr 7 2025":5*5*150/25,
    "Apr 25 2025":5*5*155/25,
    "May 2 2025":5*5*160/25,
    "May 6 2025":5*5*165/25,
    "May 29 2025":5*5*165/25,
    "Jul 2 2025":5*5*155/25,
    "Jul 17 2025":5*5*155/25,
    "Jul 25 2025":5*5*155/25,
    "Aug 20 2025":5*5*135/25,
    "Nov 19 2025":5*5*135/25,
    "Dec 3 2025":5*5*140/25,
    "Dec 10 2025":5*5*145/25,
    "Dec 17 2025":5*5*150/25,
    "Dec 28 2025":5*5*155/25,
    "Jan 14 2026":5*5*160/25,
    "Jan 22 2026":5*5*165/25,
    "Jan 28 2026":5*5*170/25,
    "Feb 14 2026":5*5*170/25,
    "Feb 18 2026":5*5*175/25,
    }

dips={
    "Jan 30 2024":9+5+5,
    # "Feb 2 2024":(6),
    # "Feb 8 2024":10+10,
    "Feb 13 2024":10+8+5,
    "Feb 16 2024":10+8+6,
    "Feb 22 2024":10+9+6,
    "Feb 27 2024":13+7+5,
    "Mar 1 2024":11+8+8,
    "Mar 6 2024":16+10+8,
    "Mar 17 2024":13+12+9,
    "Mar 27 2024":13+8+6,
    "Apr 11 2024":15+10+8,
    # "Apr 17 2024":18,
    "Apr 26 2024":18+10+8,
    "Aug 23 2024":10+7+7,
    "Sep 3 2024":10+8+6,
    "Sep 19 2024":15+9+7,
    "Feb 24 2025":8+8+5,
    "Feb 28 2025":8+6+4,
    "Mar 6 2025":8+7+5,
    "Mar 11 2025":10+7+5,
    "Mar 21 2025":10+7+5,
    "Mar 24 2025":10+7+5,
    "Mar 31 2025":10+8+6,
    "Apr 4 2025":12+10+7,
    "Apr 12 2025":12+8+5,
    "Apr 23 2025":15+10+5,
    "Apr 26 2025":12+8+8,
    "Apr 30 2025":15+9+8,
    "May 4 2025":13+9+8,
    "May 7 2025":10+10+10,
    "May 14 2025":16+10+9,
    "May 28 2025":16+10+8,
    "Jun 12 2025":16+9+5,
    "Jun 18 2025":17+12+10,
    "Jun 24 2025":17+13+10,
    "Jun 30 2025":19+13+11,
    "Jul 8 2025":20+12+11,
    "Jul 12 2025":16+10+7,
    "Jul 23 2025":16+10+7,
    "Jul 29 2025":17+13+8,
    "Aug 19 2025":16+8+6,
    "Sep 16 2025":17+10+8,
    "Oct 8 2025":16+12+9,
    "Nov 6 2025":17+10+6,
    "Nov 17 2025":12+9+7,
    "Nov 21 2025":17+11+8,
    "Dec 4 2025":20+11+10,
    "Dec 12 2025":8+8+8,
    "Dec 20 2025":15+8+7,
    "Dec 29 2025":9+9+5,
    "Jan 9 2026":12+10+7,
    }

pullups={
    "Jan 18 2023":5+5+3,
    "Apr 30 2023":(6),
    "Jul 24 2023":5+4,
    "Jul 27 2023":5+4,
    "Jan 24 2024":(3),
    "Jan 26 2024":4,
    "Jan 31 2024":3+3,
    "Feb 2 2024":5,
    "Feb 8 2024":5+4,
    "Feb 14 2024":6+4+2,
    "Feb 20 2024":5+5+4,
    "Feb 23 2024":(3),
    "Feb 28 2024":3+3+2,
    "Feb 29 2024":6+5+5+4,
    "Mar 4 2024":3+3+2,
    "Mar 7 2024":3+3+3+3+3+2+2,
    "Mar 19 2024":6+4+4+4,
    "Mar 28 2024":5+4+4+3+2,
    "Apr 12 2024":5+4+4,
    "Aug 20 2024":4+4+4,
    "Aug 28 2024":4+4+2,
    "Sep 24 2024":5+3+3,
    "Oct 2 2024":4+3+2+2+2,
    "Feb 25 2025":3+3+2,
    "Mar 4 2025":4+4+3,
    "Mar 20 2025":3+3+3,
    "Mar 29 2025":7+4+3,
    "Apr 3 2025":6+4+3,
    "Apr 11 2025":6+4+3,
    "Apr 25 2025":7+4+3,
    "Apr 28 2025":7+4+3,
    "May 2 2025":8+4+3,
    "May 6 2025":9+4+4,
    "May 9 2025":9+5+3,
    "May 29 2025":9+6+4,
    "Jun 27 2025":7+5+5,
    "Jul 2 2025":8+5+4,
    "Jul 10 2025":6+4+2,
    "Jul 17 2025":6+5+4,
    "Jul 25 2025":8+5+4,
    "Aug 5 2025":9+5+5,
    "Aug 20 2025":5+4+4,
    "Aug 29 2025":10+5+2,
    "Sep 3 2025":10+6+4,
    "Sep 17 2025":6+6+4,
    "Oct 8 2025":12+6+3,
    "Nov 8 2025":3+3+3,
    "Nov 19 2025":6+5+3,
    "Dec 3 2025":5+4+4,
    "Dec 10 2025":6+3+3,
    "Dec 17 2025":3+4+3,
    "Dec 28 2025":4+4+3,
    "Jan 8 2026":4+4+3,
    "Jan 14 2026":5+4+4,
    }

chinups={
    # "Feb 2 2024":4+3,
    "Feb 9 2024":7+4+4,
    "Feb 15 2024":8+4+3,
    "Feb 21 2024":8+5+3,
    "Feb 26 2024":11+7+6,
    "Mar 5 2024":12+7+5,
    "Mar 8 2024":11+5+5,
    "Mar 20 2024":11+6+5,
    # "Apr 5 2024":18, #didn't track individual sets
    "Apr 10 2024":11+6+5,
    "Apr 16 2024":12+7+4,
    "Apr 22 2024":12+7+4,
    "Aug 22 2024":8+5+4,
    "Aug 29 2024":8+6+4,
    "Sep 18 2024":10+6+4,
    "Oct 7 2024":8+5+4,
    "Mar 10 2025":9+5+3,
    "Mar 18 2025":9+5+3,
    "Apr 7 2025":8+4+3,
    "Jul 30 2025":5+4+4,
    "Nov 18 2025":8+5+4,
    }

bodyweight={
    "Jul 27 2016":155,
    "Aug 26 2016":145,
    "Dec 18 2016":139,
    "Feb 28 2017":157,
    "Oct 23 2017":166,
    "Dec 4 2017":166,
    "Jul 8 2018":150,
    "Dec 3 2018":158,
    "Apr 1 2020":175, #approx (ballpark estimate from memory)
    "Sep 20 2020":166,
    "Jan 21 2021":150,
    "Mar 28 2024":170.8,
    "Apr 16 2024":168.8,
    # "Apr 22 2024":167.7,
    "Apr 26 2024":171.6,
    "May 29 2024":168.2,
    "Oct 17 2024":168.4,
    "Feb 11 2025":166.5,
    "Mar 24 2025":161.2,
    "Mar 31 2025":157.1,
    # "Apr 3 2025":160, #was a full stomach
    # "Apr 8 2025":157.4,
    # "Apr 11 2025":156.4,
    # "Apr 25 2025":158.8,
    # "Apr 28 2025":156.9,
    # "May 2 2025":158.8, #after big breakfast
    # "May 6 2025":157.5,
    # # "Jun 2 2025":156,
    # "Jun 4 2025":158,
    # "Jun 12 2025":156,
    "Jun 18 2025":157,
    "Jul 8 2025":153.5,
    # "Jul 17 2025":154.5,
    # "Jul 23 2025":154,
    # "Jul 25 2025":154.5,
    "Jul 29 2025":152.5,
    # "Aug 5 2025":152.5,
    # "Aug 27 2025":152, #full stomach
    # "Aug 29 2025":150.5,
    "Sep 16 2025":150.5,
    "Oct 2 2025":147.5,
    "Nov 6 2025":150,
    # "Nov 19 2025":149,
    # "Nov 21 2025":150,
    "Dec 3 2025":149,
    # "Dec 8 2025":151,
    "Dec 12 2025":155,
    # "Dec 15 2025":155,
    # "Dec 17 2025":154,
    "Dec 21 2025":155,
    "Jan 6 2026":159,
    # "Jan 7 2026":158,
    # "Jan 9 2026":160,
    # "Jan 12 2026":158,
    # "Jan 14 2026":161,
    # "Jan 16 2026":160,
    # "Jan 20 2026":159,
    "Jan 22 2026":161,
    "Jan 31 2026":158,
    }

legextension15x3={
    "Dec 8 2025": 15*3*110/45,
    "Dec 15 2025": 15*3*125/45,
    "Jan 6 2026": 15*3*130/45,
    "Jan 19 2026": 15*3*135/45,
    "Feb 2 2026": 15*3*140/45,
    "Feb 10 2026": 15*3*150/45,
    "Feb 16 2026": 15*3*155/45,
    }

hamstringcurl15x3={
    "Dec 8 2025": 15*3*90/45,
    "Dec 15 2025": 15*3*100/45,
    "Jan 12 2026": 15*3*115/45,
    }

# lateralraise15x3={
#     "Dec 8 2025": 15*3*15/45, #fake number, but approx equivalent strength
#     "Dec 25 2025": 15*3*20/45,
#     "Jan 31 2026": 15*3*24/45,
#     }

cablecrunch15x3={
    "Dec 16 2025": 15*3*85/45,
    "Jan 7 2026": 15*3*90/45,
    "Jan 11 2026": 15*3*95/45,
    "Jan 15 2026": 15*3*100/45,
    "Jan 29 2026": 15*3*105/45,
    "Jan 31 2026": 15*3*110/45,
    "Feb 1 2026": 15*3*115/45,
    "Feb 3 2026": 15*3*120/45,
    }

shrugs10x3={
    "Dec 9 2025": 10*3*120/30,
    "Dec 16 2025": 10*3*130/30,
    "Dec 22 2025": 10*3*150/30,
    "Jan 20 2026": 10*3*160/30,
    }
rowmachine10x3={
    "Dec 3 2025": 10*3*110/30,
    "Dec 10 2025": 10*3*115/30,
    "Dec 17 2025": 10*3*120/30,
    "Jan 14 2026": 10*3*125/30,
    "Jan 22 2026": 10*3*130/30,
    "Jan 28 2026": 10*3*135/30,
    "Jan 28 2026": 10*3*135/30,
    "Feb 18 2026": 10*3*145/30,
    }

xrange=[3+365*9,170+365*11]
colours=['#2ca02c','#d62728','#1f77b4','#9467bd','#ff7f0e','#7f7f7f','#9e9ac8','#17becf','#8c564b']

ax[0].plot(date_nums(squats.keys()),squats.values(),label='squat',color=colours[0])
ax[0].plot(date_nums(deadlifts.keys()),deadlifts.values(),label='deadlifts',color=colours[1])
ax[0].axvline(x=365, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*2, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*3, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*4, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*5, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*6, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*7, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*8, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*9, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*10, linestyle='--', linewidth=1, color='grey')
ax[0].axvline(x=365*11, linestyle='--', linewidth=1, color='grey')
ax[0].set_ylabel("[lb]")
ax[0].legend(loc="upper left")
if 'xrange' in locals() or 'xrange' in globals():
    ax[0].set_xlim(xrange)

ax[1].plot(date_nums(bench.keys()),bench.values(),label='bench',color=colours[2])
ax[1].plot(date_nums(row.keys()),row.values(),label='row',color=colours[3])
ax[1].axvline(x=365, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*2, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*3, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*4, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*5, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*6, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*7, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*8, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*9, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*10, linestyle='--', linewidth=1, color='grey')
ax[1].axvline(x=365*11, linestyle='--', linewidth=1, color='grey')
ax[1].set_ylabel("[lb]")
ax[1].legend(loc="upper left")
if 'xrange' in locals() or 'xrange' in globals():
    ax[1].set_xlim(xrange)

ax[2].plot(date_nums(ohp.keys()),ohp.values(),label='ohp',color=colours[4])
ax[2].axvline(x=365, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*2, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*3, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*4, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*5, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*6, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*7, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*8, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*9, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*10, linestyle='--', linewidth=1, color='grey')
ax[2].axvline(x=365*11, linestyle='--', linewidth=1, color='grey')
ax[2].set_ylabel("[lb]")
ax[2].legend(loc="upper left")
if 'xrange' in locals() or 'xrange' in globals():
    ax[2].set_xlim(xrange)



ax[3].plot(date_nums(dips.keys()),dips.values(),label='dips',color=colours[5])
ax[3].plot(date_nums(pullups.keys()),pullups.values(),label='pullups',color=colours[6])
ax[3].plot(date_nums(chinups.keys()),chinups.values(),label='chinups',color=colours[7])
ax[3].axvline(x=365, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*2, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*3, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*4, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*5, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*6, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*7, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*8, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*9, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*10, linestyle='--', linewidth=1, color='grey')
ax[3].axvline(x=365*11, linestyle='--', linewidth=1, color='grey')
ax[3].set_ylabel("[reps]")
ax[3].legend(loc="upper left")
if 'xrange' in locals() or 'xrange' in globals():
    ax[3].set_xlim(xrange)

ax[4].plot(date_nums(bodyweight.keys()),bodyweight.values(),label='body weight',color=colours[8])
ax[4].axvline(x=365, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*2, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*3, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*4, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*5, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*6, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*7, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*8, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*9, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*10, linestyle='--', linewidth=1, color='grey')
ax[4].axvline(x=365*11, linestyle='--', linewidth=1, color='grey')
ax[4].set_xlabel("Days since Jan 2015")
ax[4].set_ylabel("[lb]")
ax[4].legend(loc="lower right")
if 'xrange' in locals() or 'xrange' in globals():
    ax[4].set_xlim(xrange)
plt.savefig("lifting2.png")
plt.show()


######################################

#PLOTTING ONLY WORKOUTS THAT SUCCESSFULLY REACHED 5x5
f, ax = plt.subplots(1, 1, figsize=(12, 9))
plt.subplots_adjust(hspace=0.7)
xrange=[19+365*9,70+365*11]
# colours=['#2ca02c','#d62728','#1f77b4','#9467bd','#ff7f0e','#7f7f7f','#9e9ac8','#17becf','#8c564b']
colours = [
    '#2ca02c', '#d62728', '#1f77b4', '#9467bd', '#ff7f0e',
    '#7f7f7f', '#9e9ac8', '#17becf', '#8c564b',
    '#e377c2', '#c49c94', '#aec7e8', '#ff9896', '#98df8a',
    '#bcbd22'
]
ax.plot(date_nums(squats5x5.keys()),squats5x5.values(),label='squat 5x5',color=colours[0])
ax.plot(date_nums(deadlifts5x5.keys()),deadlifts5x5.values(),label='deadlift 3x5',color=colours[1])
ax.plot(date_nums(bench5x5.keys()),bench5x5.values(),label='bench 5x5',color=colours[2])
ax.plot(date_nums(row5x5.keys()),row5x5.values(),label='row 5x5',color=colours[3])
ax.plot(date_nums(ohp5x5.keys()),ohp5x5.values(),label='ohp 5x5',color=colours[4])
ax.plot(date_nums(legextension15x3.keys()),legextension15x3.values(),label='leg ext. 3x15',color=colours[6])
ax.plot(date_nums(hamstringcurl15x3.keys()),hamstringcurl15x3.values(),label='ham. curl 3x15',color=colours[7])
# ax.plot(date_nums(lateralraise15x3.keys()),lateralraise15x3.values(),label='lateral raise',color=colours[8])
ax.plot(date_nums(cablecrunch15x3.keys()),cablecrunch15x3.values(),label='crunch 3x15',color=colours[9])
ax.plot(date_nums(shrugs10x3.keys()),shrugs10x3.values(),label='shrug 3x10',color=colours[10])
ax.plot(date_nums(rowmachine10x3.keys()),rowmachine10x3.values(),label='row mach. 3x10',color=colours[11])
ax.plot(date_nums(bodyweight.keys()),bodyweight.values(),label='bw',color=colours[5])
ax.axvline(x=365, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*2, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*3, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*4, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*5, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*6, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*7, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*8, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*9, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*10, linestyle='--', linewidth=1, color='grey')
ax.axvline(x=365*11, linestyle='--', linewidth=1, color='grey')
ax.set_ylabel("[lb]")
ax.set_xlabel("Days since Jan 1, 2015")
# ax.legend(loc="lower right")
ax.legend(loc="upper center",frameon=False,ncol=6)
if 'xrange' in locals() or 'xrange' in globals():
    ax.set_xlim(xrange)
ax.set_title("Weight lifted for reps")
ax.grid(True)
for x, y in zip(date_nums(squats5x5.keys()), squats5x5.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,-10), ha='center', fontsize=6,color=colours[0])
for x, y in zip(date_nums(deadlifts5x5.keys()), deadlifts5x5.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[1])
for x, y in zip(date_nums(bench5x5.keys()), bench5x5.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[2])
for x, y in zip(date_nums(row5x5.keys()), row5x5.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,-10), ha='center', fontsize=6,color=colours[3])
for x, y in zip(date_nums(ohp5x5.keys()), ohp5x5.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[4])
for x, y in zip(date_nums(legextension15x3.keys()), legextension15x3.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[6])
for x, y in zip(date_nums(hamstringcurl15x3.keys()), hamstringcurl15x3.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[7])
# for x, y in zip(date_nums(lateralraise15x3.keys()), lateralraise15x3.values()):
#     ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[8])
for x, y in zip(date_nums(cablecrunch15x3.keys()), cablecrunch15x3.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[9])
for x, y in zip(date_nums(shrugs10x3.keys()), shrugs10x3.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[10])
for x, y in zip(date_nums(rowmachine10x3.keys()), rowmachine10x3.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[11])
for x, y in zip(date_nums(bodyweight.keys()), bodyweight.values()):
    ax.annotate(int(y), (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=6,color=colours[5])
plt.savefig("lifting.png")