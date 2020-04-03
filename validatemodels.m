    function [FITje_val,dirdmd_val,x]=validatemodels(sys_red, Inputs_val,Outputs_val,r,dirdmd_val,f)
    
    FITje_val=zeros(2, length(sys_red));
    OMEGA_val={};
    DAMPING_val={};
    x=cell(r,1);
    
    if ~exist(dirdmd_val,'dir') 
            mkdir(dirdmd_val);
    end
    
    for si=1:length(sys_red)
        part=3; subpart=1; [f]= MPC_progress(part,subpart,f,si,r);
        warning off 
        [FITje_val,OMEGA,DAMPING,fig1,x]=evaluatemodel(sys_red,si,Inputs_val,Outputs_val,FITje_val,OMEGA_val,DAMPING_val,'validation',x);
        export_fig(fig1,strcat(dirdmd_val,'/image',num2str(20000+si)),'-nocrop','-m2')
        warning on
        close all
    end
        
    if r==length(sys_red)

        [fig200]=VAFpermodes(FITje_val,r,{});
        export_fig(fig200,strcat(dirdmd_val,'/image',num2str(20000+length(sys_red)+1)),'-nocrop','-m2')
        close all
    
    else
        
        Xd=[1;2;3;4];
        [fig200]=VAFpermodes(FITje_val,r,{});
        export_fig(fig200,strcat(dirdmd_val,'/image',num2str(20000+length(sys_red)+1)),'-nocrop','-m2')
        close all
    
    end
end
    