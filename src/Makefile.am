noinst_LIBRARIES = libftag.a libdblf84.a
libftag_a_SOURCES = dnasubmodel.cpp dnasubmodel.h matrix.h matrix_impl.h ctmarkov.cpp ctmarkov.h contextmodel.cpp contextmodel.h contextmodeltable.cpp contextmodeltable.h ftagmodel.cpp ftagmodel_forward.cpp ftagmodel_viterbi.cpp ftagmodel_backward.cpp ftagmodel.h ftagtable_impl.h ftagtable.h ftagtablestructs.h ftagstates.h transitionmodel.cpp transitionmodel.h ftagparams.cpp ftagparams.h ftagmodel_bw.cpp em.h transitionestimator.h transitionestimator.cpp emissionestimator.h emissionestimator.cpp ftagtrain.cpp ftagtrain.h ftaggen.h ftaggen.cpp outputqueue.h outputqueue.cpp agdoptimizer.cpp agdoptimizer.h ftagargs.cpp ftagargs.h axtpair.cpp axtpair.h axtreader.cpp axtreader.h pairs.cpp pairs.h 

libdblf84_a_SOURCES = dblf84/dblf84.cpp dblf84/dblf84.h dblf84/dblf84_dec.h dblf84/dblf84_111.cpp dblf84/dblf84_112.cpp dblf84/dblf84_113.cpp dblf84/dblf84_114.cpp dblf84/dblf84_121.cpp dblf84/dblf84_122.cpp dblf84/dblf84_123.cpp dblf84/dblf84_124.cpp dblf84/dblf84_131.cpp dblf84/dblf84_132.cpp dblf84/dblf84_133.cpp dblf84/dblf84_134.cpp dblf84/dblf84_141.cpp dblf84/dblf84_142.cpp dblf84/dblf84_143.cpp dblf84/dblf84_144.cpp dblf84/dblf84_211.cpp dblf84/dblf84_212.cpp dblf84/dblf84_213.cpp dblf84/dblf84_214.cpp dblf84/dblf84_221.cpp dblf84/dblf84_222.cpp dblf84/dblf84_223.cpp dblf84/dblf84_224.cpp dblf84/dblf84_231.cpp dblf84/dblf84_232.cpp dblf84/dblf84_233.cpp dblf84/dblf84_234.cpp dblf84/dblf84_241.cpp dblf84/dblf84_242.cpp dblf84/dblf84_243.cpp dblf84/dblf84_244.cpp dblf84/dblf84_311.cpp dblf84/dblf84_312.cpp dblf84/dblf84_313.cpp dblf84/dblf84_314.cpp dblf84/dblf84_321.cpp dblf84/dblf84_322.cpp dblf84/dblf84_323.cpp dblf84/dblf84_324.cpp dblf84/dblf84_331.cpp dblf84/dblf84_332.cpp dblf84/dblf84_333.cpp dblf84/dblf84_334.cpp dblf84/dblf84_341.cpp dblf84/dblf84_342.cpp dblf84/dblf84_343.cpp dblf84/dblf84_344.cpp dblf84/dblf84_411.cpp dblf84/dblf84_412.cpp dblf84/dblf84_413.cpp dblf84/dblf84_414.cpp dblf84/dblf84_421.cpp dblf84/dblf84_422.cpp dblf84/dblf84_423.cpp dblf84/dblf84_424.cpp dblf84/dblf84_431.cpp dblf84/dblf84_432.cpp dblf84/dblf84_433.cpp dblf84/dblf84_434.cpp dblf84/dblf84_441.cpp dblf84/dblf84_442.cpp dblf84/dblf84_443.cpp dblf84/dblf84_444.cpp

bin_PROGRAMS = smtest cmtest tmtest ftagtest tabletest genpair traintest opttest axtscan simulate axtsample train align homolog compare

smtest_SOURCES = dnasubmodeltest.cpp 
smtest_LDADD = libftag.a libdblf84.a

cmtest_SOURCES = contextmodeltest.cpp
cmtest_LDADD = libftag.a libdblf84.a

tmtest_SOURCES = transitionmodeltest.cpp
tmtest_LDADD = libftag.a libdblf84.a

ftagtest_SOURCES = ftagmodeltest.cpp 
ftagtest_LDADD = libftag.a libdblf84.a

tabletest_SOURCES = ftagtabletest.cpp
tabletest_LDADD = libftag.a libdblf84.a

genpair_SOURCES = genpair.cpp
genpair_LDADD = libftag.a libdblf84.a

traintest_SOURCES = ftagtraintest.cpp
traintest_LDADD = libftag.a libdblf84.a

opttest_SOURCES = optimizetest.cpp
opttest_LDADD = libftag.a libdblf84.a

axtscan_SOURCES = axtscan.cpp
axtscan_LDADD = libftag.a libdblf84.a

simulate_SOURCES = simulate.cpp
simulate_LDADD = libftag.a libdblf84.a

axtsample_SOURCES = axtsample.cpp
axtsample_LDADD = libftag.a libdblf84.a

train_SOURCES = train.cpp
train_LDADD = libftag.a libdblf84.a

align_SOURCES = align.cpp
align_LDADD = libftag.a libdblf84.a

homolog_SOURCES = homolog.cpp
homolog_LDADD = libftag.a libdblf84.a

compare_SOURCES = compare.cpp
compare_LDADD = libftag.a libdblf84.a