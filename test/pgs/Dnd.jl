D2ds = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("C_2^1", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("S_4^1", [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]),
    Symel("S_4^3", [0.0 1.0 0.0; -1.0 0.0 -0.0; 0.0 0.0 -1.0]),
    Symel("sigmad_2", [0.0 -1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [0.0 1.0 0.0; 1.0 0.0 -0.0; 0.0 -0.0 1.0]),
    Symel("C_2'(2)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(1)", [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0])]

D3ds = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("i", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_3^1", [-0.4999999999999998 -0.8660254037844387 0.0; 0.8660254037844387 -0.4999999999999998 0.0; 0.0 0.0 1.0]),
    Symel("C_3^2", [-0.5000000000000003 0.8660254037844384 0.0; -0.8660254037844384 -0.5000000000000003 0.0; 0.0 0.0 1.0]),
    Symel("S_6^1", [0.5000000000000001 -0.8660254037844386 0.0; 0.8660254037844386 0.5000000000000001 0.0; 0.0 0.0 -1.0]),
    Symel("S_6^5", [0.49999999999999944 0.8660254037844392 0.0; -0.8660254037844392 0.49999999999999944 0.0; 0.0 0.0 -1.0]),
    Symel("sigmad_2", [-1.0 0.0 -0.0; 0.0 1.0 0.0; -0.0 0.0 1.0]),
    Symel("sigmad_3", [0.49999999999999933 -0.8660254037844392 0.0; -0.8660254037844392 -0.4999999999999998 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [0.5000000000000008 0.8660254037844383 0.0; 0.8660254037844383 -0.5000000000000009 -0.0; 0.0 -0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-0.5000000000000002 -0.8660254037844387 0.0; -0.8660254037844387 0.5000000000000009 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(3)", [-0.4999999999999991 0.8660254037844394 0.0; 0.8660254037844394 0.4999999999999998 0.0; 0.0 0.0 -1.0])]

D4ds = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("C_4^1", [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0]),
    Symel("C_2^1", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("C_4^3", [0.0 1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0]),
    Symel("S_8^1", [0.7071067811865476 -0.7071067811865475 0.0; 0.7071067811865475 0.7071067811865476 0.0; 0.0 0.0 -1.0]),
    Symel("S_8^3", [-0.7071067811865474 -0.7071067811865477 -0.0; 0.7071067811865477 -0.7071067811865474 0.0; 0.0 0.0 -1.0]),
    Symel("S_8^5", [-0.7071067811865479 0.7071067811865471 0.0; -0.7071067811865471 -0.7071067811865479 -0.0; 0.0 0.0 -1.0]),
    Symel("S_8^7", [0.707106781186547 0.707106781186548 0.0; -0.707106781186548 0.707106781186547 0.0; 0.0 0.0 -1.0]),
    Symel("sigmad_2", [-0.7071067811865475 0.7071067811865477 -0.0; 0.7071067811865477 0.7071067811865475 0.0; -0.0 0.0 1.0]),
    Symel("sigmad_3", [-0.7071067811865479 -0.7071067811865472 0.0; -0.7071067811865472 0.7071067811865479 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_4", [0.707106781186547 -0.7071067811865481 0.0; -0.7071067811865481 -0.7071067811865479 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [0.7071067811865481 0.707106781186547 0.0; 0.707106781186547 -0.7071067811865483 -0.0; 0.0 -0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(1)", [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(2)", [0.0 -1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 -1.0])]

D5ds = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("i", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_5^1", [0.30901699437494745 -0.9510565162951535 0.0; 0.9510565162951535 0.30901699437494745 0.0; 0.0 0.0 1.0]),
    Symel("C_5^2", [-0.8090169943749473 -0.5877852522924731 0.0; 0.5877852522924731 -0.8090169943749473 0.0; 0.0 0.0 1.0]),
    Symel("C_5^3", [-0.8090169943749475 0.587785252292473 0.0; -0.587785252292473 -0.8090169943749475 0.0; 0.0 0.0 1.0]),
    Symel("C_5^4", [0.30901699437494734 0.9510565162951535 0.0; -0.9510565162951535 0.30901699437494734 0.0; 0.0 0.0 1.0]), 
    Symel("S_10^1", [0.8090169943749475 -0.5877852522924731 0.0; 0.5877852522924731 0.8090169943749475 0.0; 0.0 0.0 -1.0]),
    Symel("S_10^3", [-0.3090169943749474 -0.9510565162951536 -0.0; 0.9510565162951536 -0.3090169943749474 0.0; 0.0 0.0 -1.0]),
    Symel("S_10^7", [-0.30901699437494756 0.9510565162951538 0.0; -0.9510565162951538 -0.30901699437494756 -0.0; 0.0 0.0 -1.0]),
    Symel("S_10^9", [0.8090169943749477 0.5877852522924735 0.0; -0.5877852522924735 0.8090169943749477 0.0; 0.0 0.0 -1.0]),
    Symel("sigmad_2", [-0.30901699437494745 0.9510565162951536 -0.0; 0.9510565162951536 0.30901699437494745 0.0; -0.0 0.0 1.0]),
    Symel("sigmad_3", [-1.0 -0.0 -0.0; -0.0 1.0 -0.0; -0.0 -0.0 1.0]),
    Symel("sigmad_4", [-0.30901699437494745 -0.9510565162951536 0.0; -0.9510565162951536 0.30901699437494745 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_5", [0.8090169943749473 -0.5877852522924732 0.0; -0.5877852522924732 -0.8090169943749472 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [0.8090169943749475 0.587785252292473 0.0; 0.587785252292473 -0.8090169943749472 -0.0; 0.0 -0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-0.8090169943749473 0.5877852522924734 0.0; 0.5877852522924734 0.8090169943749481 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(3)", [0.3090169943749479 -0.951056516295154 0.0; -0.951056516295154 -0.3090169943749471 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(4)", [0.30901699437494834 0.9510565162951539 0.0; 0.9510565162951539 -0.30901699437494745 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(5)", [-0.8090169943749475 -0.5877852522924731 0.0; -0.5877852522924731 0.8090169943749481 0.0; 0.0 0.0 -1.0])]

D6ds = [Symel("E", [1 0 0; 0 1 0; 0 0 1]),
    Symel("C_6^1", [0.5000000000000001 -0.8660254037844386 0.0; 0.8660254037844386 0.5000000000000001 0.0; 0.0 0.0 1.0]),
    Symel("C_3^1", [-0.4999999999999998 -0.8660254037844388 0.0; 0.8660254037844388 -0.4999999999999998 0.0; 0.0 0.0 1.0]),
    Symel("C_2^1", [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0]),
    Symel("C_3^2", [-0.5000000000000006 0.8660254037844385 0.0; -0.8660254037844385 -0.5000000000000006 0.0; 0.0 0.0 1.0]),
    Symel("C_6^5", [0.49999999999999944 0.8660254037844392 0.0; -0.8660254037844392 0.49999999999999944 0.0; 0.0 0.0 1.0]),
    Symel("S_12^1", [0.8660254037844387 -0.49999999999999994 0.0; 0.49999999999999994 0.8660254037844387 0.0; 0.0 0.0 -1.0]),
    Symel("S_4^1", [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0]),
    Symel("S_12^5", [-0.8660254037844385 -0.5000000000000006 -0.0; 0.5000000000000006 -0.8660254037844385 0.0; 0.0 0.0 -1.0]),
    Symel("S_12^7", [-0.866025403784439 0.49999999999999944 0.0; -0.49999999999999944 -0.866025403784439 -0.0; 0.0 0.0 -1.0]),
    Symel("S_4^3", [0.0 1.0 0.0; -1.0 0.0 -0.0; 0.0 0.0 -1.0]),
    Symel("S_12^11", [0.8660254037844383 0.500000000000001 0.0; -0.500000000000001 0.8660254037844383 0.0; 0.0 0.0 -1.0]),
    Symel("sigmad_2", [2.220446049250313e-16 1.0 -0.0; 1.0 -2.220446049250313e-16 0.0; -0.0 0.0 1.0]),
    Symel("sigmad_3", [-0.8660254037844384 0.5000000000000004 -0.0; 0.5000000000000004 0.8660254037844384 0.0; -0.0 0.0 1.0]),
    Symel("sigmad_4", [-0.8660254037844393 -0.4999999999999995 0.0; -0.4999999999999995 0.8660254037844389 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_5", [-1.1102230246251565e-15 -1.0000000000000002 0.0; -1.0000000000000002 7.771561172376096e-16 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_6", [0.8660254037844379 -0.5000000000000012 0.0; -0.5000000000000012 -0.8660254037844384 0.0; 0.0 0.0 1.0]),
    Symel("sigmad_1", [0.8660254037844393 0.49999999999999883 0.0; 0.49999999999999883 -0.8660254037844393 -0.0; 0.0 -0.0 1.0]),
    Symel("C_2'(1)", [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(2)", [-0.4999999999999998 0.8660254037844388 0.0; 0.8660254037844388 0.4999999999999998 0.0; 0.0 0.0 -1.0]),
    Symel("C_2'(3)", [-0.5000000000000004 -0.8660254037844385 0.0; -0.8660254037844385 0.5000000000000007 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(1)", [0.5000000000000002 0.8660254037844386 0.0; 0.8660254037844386 -0.5000000000000001 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(2)", [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]),
    Symel("C_2''(3)", [0.4999999999999998 -0.8660254037844392 0.0; -0.8660254037844392 -0.49999999999999933 0.0; 0.0 0.0 -1.0])]

D2dcn = ["A1", "A2", "B1", "B2", "E"]
D2dct = [1.0  1.0  1.0  1.0  1.0;
         1.0  1.0  1.0 -1.0 -1.0;
         1.0 -1.0  1.0  1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0;
         2.0  0.0 -2.0  0.0  0.0]
D3dcn = ["A1g", "A2g", "Eg", "A1u", "A2u", "Eu"]
D3dct = [1.0  1.0  1.0  1.0  1.0  1.0;
         1.0  1.0 -1.0  1.0  1.0 -1.0;
         2.0 -1.0  0.0  2.0 -1.0  0.0;
         1.0  1.0  1.0 -1.0 -1.0 -1.0;
         1.0  1.0 -1.0 -1.0 -1.0  1.0;
         2.0 -1.0  0.0 -2.0  1.0  0.0]
D4dcn = ["A1", "A2", "B1", "B2", "E1", "E2", "E3"]
D4dct = [1.0  1.0  1.0  1.0  1.0  1.0  1.0;
         1.0  1.0  1.0  1.0  1.0 -1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0  1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0;
         2.0  rt2  0.0 -rt2 -2.0  0.0  0.0;
         2.0  0.0 -2.0  0.0  2.0  0.0  0.0;
         2.0 -rt2  0.0  rt2 -2.0  0.0  0.0]
D5dcn = ["A1g", "A2g", "E1g", "E2g", "A1u", "A2u", "E1u", "E2u"]
D5dct = [1.0  1.0       1.0        1.0  1.0  1.0       1.0        1.0;
         1.0  1.0       1.0       -1.0  1.0  1.0       1.0       -1.0;
         2.0  twocos72  twocos144  0.0  2.0  twocos72  twocos144  0.0;
         2.0  twocos144 twocos72   0.0  2.0  twocos144 twocos72   0.0;
         1.0  1.0       1.0        1.0 -1.0 -1.0      -1.0       -1.0;
         1.0  1.0       1.0       -1.0 -1.0 -1.0      -1.0        1.0;
         2.0  twocos72  twocos144  0.0 -2.0 -twocos72 -twocos144  0.0;
         2.0  twocos144 twocos72   0.0 -2.0 -twocos144 -twocos72   0.0]
D6dcn = ["A1", "A2", "B1", "B2", "E1", "E2", "E3", "E4", "E5"]
D6dct = [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
         1.0  1.0  1.0  1.0  1.0  1.0  1.0 -1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0  1.0 -1.0;
         1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0;
         2.0  rt3  1.0  0.0 -1.0 -rt3 -2.0  0.0  0.0;
         2.0  1.0 -1.0 -2.0 -1.0  1.0  2.0  0.0  0.0;
         2.0  0.0 -2.0  0.0  2.0  0.0 -2.0  0.0  0.0;
         2.0 -1.0 -1.0  2.0 -1.0 -1.0  2.0  0.0  0.0;
         2.0 -rt3  1.0  0.0 -1.0  rt3 -2.0  0.0  0.0]

