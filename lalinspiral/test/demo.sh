make ChirpSpace
make CoarseTest
make CoarseTest2
./CoarseTest > CoarseTest.out
./CoarseTest2 > CoarseTest2.out
./ChirpSpace > ChirpSpace.out
xmgrace -par CoarseTest.par ChirpSpace.out CoarseTest.out &
xmgrace -par CoarseTest2.par ChirpSpace.out CoarseTest2.out &
make SpaceCovering
./SpaceCovering --template EOB --grid-spacing Hexagonal ; xmgrace SpaceCovering.out ChirpSpace.out -p SpaceCovering.par
