target_sources(MarDyn
    PRIVATE
        AutoPasClass.cpp
        iteratePairwiseLJFunctorAVXNoShiftMix.cpp
        iteratePairwiseLJFunctorAVXNoShiftNoMix.cpp
        iteratePairwiseLJFunctorAVXShiftMix.cpp
        iteratePairwiseLJFunctorAVXShiftNoMix.cpp
        iteratePairwiseLJFunctorNoShiftMix.cpp
        iteratePairwiseLJFunctorNoShiftNoMix.cpp
        iteratePairwiseLJFunctorShiftMix.cpp
        iteratePairwiseLJFunctorShiftNoMix.cpp
        iteratePairwiseLJFunctorSVENoShiftMix.cpp
        iteratePairwiseLJFunctorSVENoShiftNoMix.cpp
        iteratePairwiseLJFunctorSVEShiftMix.cpp
        iteratePairwiseLJFunctorSVEShiftNoMix.cpp
    )

