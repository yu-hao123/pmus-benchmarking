yalmipmainpath = 'YALMIP-master';

addpath(genpath(yalmipmainpath))
addpath(genpath(strcat(yalmipmainpath,'\extras')))
addpath(genpath(strcat(yalmipmainpath,'\solvers')))
addpath(genpath(strcat(yalmipmainpath,'\modules')))
addpath(genpath(strcat(yalmipmainpath,'\parametric')))
addpath(genpath(strcat(yalmipmainpath,'\modules\moment')))
addpath(genpath(strcat(yalmipmainpath,'\modules\global')))
addpath(genpath(strcat(yalmipmainpath,'\modules\robust')))
addpath(genpath(strcat(yalmipmainpath,'\modules\sos')))
addpath(genpath(strcat(yalmipmainpath,'\operators')))

yalmiptest

clear yalmipmainpath