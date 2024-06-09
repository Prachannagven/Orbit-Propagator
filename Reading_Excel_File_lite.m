clc
clear

table_pull = readtable  ("Drag_Constants.xlsx", Sheet="Sheet1");

array_h0 = table2array([table_pull(:,1)]);
array_rho0 = table2array([table_pull(:,2)]);
array_H = table2array([table_pull(:,3)]);

h = 99;

h0 = lookup(h, table_pull)

rho0 = interp1(array_h0, array_rho0, h_0);

function h_0 = lookup(h, table_pull)
    for i = 1:1:36
        h_0 = table2array([table_pull(i,1)]);
        if h<h_0
            if i == 1
                h_0 = table2array([table_pull(i,1)]);
            else
                h_0 = table2array([table_pull(i-1,1)]);
            end
            break
        end
    end
end

