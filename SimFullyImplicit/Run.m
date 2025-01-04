q     = [0.1, 0.2, 0.5, 1, 2, 5, 10]';
PV = arrayfun(@(n)NumCoreSimFullyImplicit_2D(q(n), 10), 1:7);
PVdata = [q,PV];
save PVdatan PVdata