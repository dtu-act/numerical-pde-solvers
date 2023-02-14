function writeImpedanceParams(impedance_data,acc,c,boundary_type,data_filename,path_json)    
    data.comment__ = sprintf("parameters for c=%i",c);
    data.boundary_type = boundary_type;
    data.data_filename = data_filename;
    data.impedance_data = impedance_data;
    
    acc_data.phi_0_max = minmax(acc(:,1));
    acc_data.phi_0_norm = abs(round(1/acc_data.phi_0_max,2));
    acc_data.phi_1_max = minmax(acc(:,2));
    acc_data.phi_1_norm = abs(round(1/acc_data.phi_1_max,2));
    acc_data.psi0_0_max = minmax(acc(:,3));
    acc_data.psi0_0_norm = abs(round(1/acc_data.psi0_0_max,2));
    acc_data.psi1_0_max = minmax(acc(:,4));
    acc_data.psi1_0_norm = abs(round(1/acc_data.psi1_0_max,2));
    
    data.acc_data = acc_data;
        
    json = jsonencode(data,'PrettyPrint',true);
    
    fileID = fopen(path_json,'w+');
    fprintf(fileID,json);
    fclose(fileID);
end

function a = minmax(arr)
 [~,i] = max(abs(arr));
 a = arr(i); 
end