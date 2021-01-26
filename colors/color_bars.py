#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:35:44 2019

@author: tylergreen
"""

from matplotlib import colors, cm
import numpy as np

def Tyler_colors(variable):
    
    if variable == 'wind_speed': 
        
        #colors_wind = [(255,255,255),(59,245,239),(57,234,228),(54,223,218),
        #               (51,213,208),(48,201,196),(45,189,184),(40,172,168),
        #               (35,157,153),(41, 173, 13),(51, 201, 20),(117, 232, 109),
        #               (184,252,120),(240,251,131),(251,215,131),
        #               (245,186,51),(246,98,29),(238,62,27),(241,15,15),
        #               (195,14,14),(187,23,23),(162,18,18),(120,12,12),
        #               (202,83,234),(208,107,235),(211,120,236),(216,147,234),
        #               (222,170,236),(227,197,235),(229,215,233),(235,162,181),
        #               (223,145,165),(213,128,150),(191,117,136),(173,94,114),
        #               (156,75,96),(165,34,68),(158,25,59),(137,9,42)]
   
        colors_wind = [(193,247,248),(118,238,242),(97,222,228),(80,206,213),(61,176,185),(42,147,158),
                       (27,143,159),(14,141,153),(5,162,97),(9,183,54),(62,198,69),(114,214,84),(168,229,100),
                       (212,237,105),(234,229,94),(255,209,77),(255,150,56),(255,93,35),(255,37,14),(241,0,0),
                       (199,0,0),(159,0,0),(116,0,0),(117,9,64),(162,26,196),(191,55,255),(206,92,255),
                       (220,129,255),(233,159,255),(245,182,249),(252,193,238),(255,188,220),(255,171,200),(255,155,178),
                       (255,139,158),(255,124,136),(254,107,114),(244,91,99),(231,75,84),(207,59,78),(180,42,71),(166,33,69)]

        
        colors_wind = [tuple([i / 255 for i in j]) for j in colors_wind]
        
        #ticks = np.arange(4,82,2)
        ticks = [7.0,10.0,13.0,16.0,19.0,22.0,25.0,28.0,31.0,34.0,37.0,40.0,42.0,44.0,46.0,48.0,50.0,52.0,54.0,56.0,58.0,60.0,62.0,64.0,69.33,74.66,80.0,85.33,
                 90.66,96.0,100.83,105.67,110.5,115.33,120.167,125.0,130.0,135.0,140.0,145.0,150.0,155.0]
        cmap = colors.ListedColormap(colors_wind)
        cmap.set_under((1,1,1))
        cmap.set_over((153/255,24/255,66/255))
        norm = colors.BoundaryNorm(ticks,cmap.N)
        


    elif variable == 'radial_velocity':
        colors_rv =[(255,0,132),   #-64
                    (249,0,132),   #-63
                    (237,0,134),   #-62
                    (225,1,136),   #-61
                    (219,1,137),   #-60
                    (200,1,139),   #-59
                    (194,2,140),   #-58
                    (175,2,142),   #-57
                    (169,2,143),   #-56
                    (151,3,146),   #-55
                    (138,3,147),   #-54
                    (126,4,149),   #-53
                    (108,4,152),   #-52
                    (93,5,153),    #-51
                    (77,4,153),    #-50
                    (61,3,153),    #-49
                    (52,3,153),    #-48
                    (28,2,153),    #-47
                    (23,23,160),   #-46
                    (24,34,163),   #-45
                    (26,55,170),   #-44
                    (28,76,177),   #-43
                    (31,108,188),  #-42
                    (32,119,192),  #-41
                    (33,130,195),  #-40
                    (36,162,205),  #-39
                    (39,193,216),  #-38
                    (40,204,220),  #-37
                    (48,224,227),  #-36
                    (65,225,228),  #-35
                    (84,227,230),  #-34
                    (91,228,231),  #-33
                    (104,229,232), #-32
                    (117,231,234), #-31
                    (136,232,236), #-30
                    (143,234,237), #-29
                    (149,234,237), #-28
                    (169,237,240), #-27
                    (182,238,241), #-26
                    (154,242,199), #-25
                    (140,243,181), #-24
                    (100,247,126), #-23
                    (73,250,89),   #-22
                    (46,252,52),   #-21
                    (3,250,3),     #-20
                    (2,240,3),     #-19
                    (3,229,3),     #-18
                    (2,219,3),     #-17
                    (2,214,3),     #-16
                    (3,198,3),     #-15
                    (3,187,3),     #-14
                    (3,182,3),     #-13
                    (2,171,2),     #-12
                    (2,156,2),     #-11
                    (2,150,2),     #-10
                    (2,135,2),     #-9
                    (2,119,2),     #-8
                    (2,114,2),     #-7
                    (2,103,2),     #-6
                    (78,121,76),   #-5
                    (90,125,88),   #-4
                    (94,126,92),   #-3
                    (102,129,100), #-2
                    (110,132,108), #-1
                    (114,133,112), #-0
                    (138,114,129), #0
                    (138,108,122), #1
                    (136,95,108),  #2
                    (135,82,94) ,  #3
                    (133,69,79),   #4
                    (132,56,65),   #5
                    (110,0,0),     #6
                    (126,0,0),     #7
                    (132,0,1),     #8
                    (149,0,2),     #9
                    (160,0,2),     #10
                    (171,0,3),     #11
                    (176,0,3),     #12
                    (193,0,4),     #13
                    (200,2,6),     #14
                    (214,0,5),     #15
                    (227,0,6),     #16
                    (238,0,7),     #17
                    (251,60,89),   #18
                    (250,65,97),   #19
                    (250,71,105),  #20
                    (252,87,130),  #21
                    (252,98,146),  #22
                    (252,109,163), #23
                    (253,120,179), #24
                    (254,131,196), #25
                    (254,137,204), #26
                    (255,159,203), #27
                    (255,178,193), #28
                    (255,197,183), #29
                    (255,215,173), #30
                    (255,232,163), #31
                    (255,228,159), #32
                    (255,215,147), #33
                    (255,211,142), #34
                    (255,197,130), #35 
                    (255,193,125), #36
                    (255,180,113), #37
                    (255,176,108), #38
                    (255,162,96),  #39
                    (255,154,87),  #40
                    (255,138,79),  #41
                    (252,135,78),  #42
                    (241,126,72),  #43
                    (234,120,69),  #44
                    (227,114,65),  #45
                    (220,108,62),  #46
                    (213,101,58),  #47
                    (209,98,56),   #48
                    (199,89,51),   #49
                    (192,83,47),   #50
                    (185,77,44),   #51
                    (177,70,40),   #52
                    (170,64,36),   #53
                    (167,61,35),   #54
                    (156,52,29),   #55
                    (149,46,26),   #56
                    (142,40,22),   #57
                    (138,36,20),   #58
                    (128,27,15),   #59
                    (121,21,11),   #60
                    (114,15,8),    #61
                    (107,9,4),     #62
                    (0,0,0)        #63
                    ]
                   
        colors_rv = [tuple([i / 255 for i in j]) for j in colors_rv]
        
        ticks_first = np.arange(-64.5,0.5,1)
        ticks_second = np.zeros(1)
        ticks_third = np.arange(0.5,64.5,1)
        ticks = np.concatenate([ticks_first,ticks_second,ticks_third])
        cmap = colors.ListedColormap(colors_rv)
        norm = colors.BoundaryNorm(ticks,cmap.N)
   
    #Reflectivity
    elif  variable == 'z':
        colors_z = [ (105,27,26),    #59.5 - 60
                     (111,26,25),    #59 - 59.5
                     (117,25,25),    #58.5 - 59
                     (124,25,25),    #58 - 58.5
                     (131,25,25),    #57.5 - 58
                     (137,24,25),    #57 - 57.5
                     (144,23,25),    #56.5 - 57 
                     (150,24,25),    #56 - 56.5
                     (158,23,25),    #55.5 - 56
                     (165,22,26),    #55 - 55.5
                     (171,22,25),    #54.5 - 55
                     (177,21,26),    #54 - 54.5
                     (185,20,26),    #53.5 - 54
                     (192,19,27),    #53 - 53.5
                     (198,19,27),    #52.5 - 53
                     (205,18,28),    #52 - 52.5
                     (213,16,28),    #51.5 - 52
                     (220,13,29),    #51 - 51.5
                     (228,11,30),    #50.5 - 51
                     (235,12,31),    #50 - 50.5
                     (155,68,24),    #49.5 - 50
                     (159,70,24),    #49 - 49.5
                     (163,73,24),    #48.5 - 49
                     (168,76,25),    #48 - 48.5
                     (171,79,25),    #47.5 - 48
                     (175,82,26),    #47 - 47.5
                     (180,85,27),    #46.5 - 47
                     (185,88,27),    #46 - 46.5
                     (189,91,28),    #45.5 - 46
                     (193,94,29),    #45 - 45.5
                     (197,98,29),    #44.5 - 45
                     (202,101,30),   #44 - 44.5
                     (206,104,31),   #43.5 - 44
                     (211,106,31),   #43 - 43.5
                     (215,109,32),   #42.5 - 43
                     (219,112,33),   #42 - 42.5
                     (223,115,34),   #41.5 - 42
                     (228,118,35),   #41 - 41.5
                     (232,121,35),   #40.5 - 41
                     (237,125,36),   #40 - 40.5
                     (179,160,35),   #39.5 - 40
                     (182,165,37),   #39 - 39.5
                     (187,171,38),   #38.5 - 39
                     (192,176,39),   #38 - 38.5
                     (195,182,41),   #37.5 - 38
                     (199,188,42),   #37 - 37.5
                     (204,194,44),   #36.5 - 37
                     (208,199,45),   #36 - 36.5
                     (213,205,47),   #35.5 - 36
                     (218,212,48),   #35 - 35.5
                     (221,218,49),   #34.5 - 35
                     (226,224,51),   #34 - 34.5
                     (231,230,52),   #33.5 - 34
                     (234,236,53),   #33 - 33.5
                     (239,242,54),   #32.5 - 33
                     (244,246,56),   #32 - 32.5
                     (216,229,51),   #31.5 - 32
                     (188,212,46),   #31 - 31.5
                     (162,194,41),   #30.5 - 31
                     (136,175,36),   #30 - 30.5
                     (112,159,31),   #29.5 - 30
                     (89,143,27),    #29 - 29.5
                     (68,126,24),    #28.5 - 29
                     (49,110,21),    #28 - 28.5
                     (32,95,19),     #27.5 - 28
                     (20,81,17),     #27 - 27.5
                     (21,87,18),     #26.5 - 27
                     (23,92,21),     #26 - 26.5
                     (25,99,25),     #25.5 - 26
                     (27,105,28),    #25 - 25.5
                     (29,111,32),    #24.5 - 25
                     (31,118,36),    #24 - 24.5
                     (34,125,40),    #23.5 - 24
                     (36,131,44),    #23 - 23.5
                     (38,138,48),    #22.5 - 23
                     (41,145,53),    #22 - 22.5
                     (43,151,57),    #21.5 - 22
                     (46,158,61),    #21 - 21.5
                     (49,165,65),    #20.5 - 21
                     (55,155,90),    #20 - 20.5
                     (65,146,121),   #19.5 - 20
                     (76,136,152),   #19 - 19.5
                     (74,132,150),   #18.5 - 19
                     (72,128,147),   #18 - 18.5
                     (70,124,145),   #17.5 -18
                     (68,120,143),   #17 - 17.5
                     (66,115,141),   #16.5 - 17
                     (64,111,137),   #16 - 16.5
                     (62,107,135),   #15.5 - 16
                     (61,103,133),   #15 - 15.5
                     (59,99,131),    #14.5 - 15
                     (57,95,128),    #14 - 14.5
                     (56,91,125),    #13.5 - 14
                     (54,88,124),    #13 - 13.5
                     (52,83,121),    #12.5 - 13
                     (51,80,120),    #12 - 12.5
                     (49,76,118),    #11.5 - 12
                     (47,72,114),    #11 - 11.5
                     (46,68,112),    #10.5 - 11
                     (44,65,110),    #10 - 10.5
                     (42,62,108),    #9.5 - 10
                     (41,59,106),    #9 - 9.5
                     (44,61,106),    #8.5 - 9
                     (48,65,107),    #8 - 8.5
                     (53,68,108),    #7.5 - 8
                     (56,71,108),    #7 - 7.5
                     (59,76,111),    #6.5 - 7
                     (64,80,112),    #6 - 6.5
                     (68,84,113),    #5.5 - 6
                     (73,88,114),    #5 - 5.5
                     (76,91,115),    #4.5 - 5
                     (81,94,116),    #4 - 4.5
                     (86,98,117),    #3.5 - 4
                     (91,103,118),   #3 - 3.5
                     (94,106,120),   #2.5 - 3
                     (99,111,121),   #2 - 2.5
                     (104,115,122),  #1.5 - 2
                     (109,119,123),  #1 - 1.5
                     (115,124,124),  #0.5 - 1
                     (118,126,125),  #0 - 0.5
                     (123,130,126),  #-0.5 - 0
                     (128,134,127),  #-1 - -0.5
                     (133,129,128),  #-1.5 - -1
                     (138,143,129),  #-2 - -1.5
                     (143,147,131)   #-2.5 - -2
                     ]
        colors_z.reverse()
        colors_z = [tuple([i / 255 for i in j]) for j in colors_z]    
        ticks = np.arange(-2.5,60.5,0.5)
        cmap = colors.ListedColormap(colors_z)
        norm = colors.BoundaryNorm(ticks,cmap.N)

    return cmap, norm, ticks 
