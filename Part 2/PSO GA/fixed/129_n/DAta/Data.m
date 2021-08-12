function [dist_ij, node_stats]  = Get_Graph_Data()

dist_ij = [0   40  60  75  90  200 100 160 80 ;
           40  0   65  40  100 50  75  110 100;
           60  65  0   75  100 100 75  75  75 ;
           75  40  75  0   100 50  90  90  150;
           90  100 100 100 0   100 75  75  100;
           200 50  100 50  100 0   70  90  75 ;
           100 75  75  90  75  70  0   70  100;
           160 110 75  90  75  90  70  0   100;
           80  100 75  150 100 75  100 100 0  ;
           ];
    
node_stats = [0   2   1.5 4.5 3   1.5 4   2.5 3  ; % (g_i) Amount of cargo to Service
              0   1   2   1   3   2   2.5 3   0.8; % (t_i) Amount of time for Servicing
              0   1   4   1   4   3   2   5   1.5; % (ET_i)Earliest possible time
              9   4   6   2   7   5.5 5   8   4  ; % (LT_i)Latest possible time
             ];
end