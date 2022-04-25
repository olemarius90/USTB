function flowField = Phantom_small2Dtube( setup ) 

%% small 2D tube phantom, run and check signal integrity in the middle of the tube
p.btf = 60;
p.npoints = 10;
p.flowlength = 0.005;
p.tubedepth = 0.015; %0.03;
p.depthstep = 0.00015; %lambda/2 for 5 MHz
p.noFlowLines = 1; %odd number
p.max_vel = 1; % [m/s]


fields = fieldnames(setup);
for k=1:size(fields,1)
    if(isfield(p,fields{k}))
        p.(fields{k}) = setup.(fields{k});
    else
        disp([ fields{k} ' is not a valid parameter for this phantom type']);
    end
end


depthtab = (-(p.noFlowLines-1)/2:1:(p.noFlowLines-1)/2)*p.depthstep+p.tubedepth;
time_max = p.flowlength/p.max_vel;


for kk = 1:p.noFlowLines,
    currtubedepth = depthtab(kk);
    flowField(kk).timetab = linspace(0, time_max, p.npoints);
    flowField(kk).postab = p.max_vel*(flowField(kk).timetab-time_max/2).*[sind(p.btf); 0; cosd(p.btf)]+[0; 0; currtubedepth];
    flowField(kk).timetab = flowField(kk).timetab.'; 
    flowField(kk).postab = flowField(kk).postab.';
end