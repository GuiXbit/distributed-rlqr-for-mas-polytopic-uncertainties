function DELTA =CalcDelta(s,l,opts,k)
    if nargin < 3 || isempty(opts)
        opts = "default";
    end
    
    switch lower(string(opts))
        case "default"
            DELTA = (2*rand(s,l) - 1)*0.5;
        case "sin"
            DELTA = sin(k*pi/4) * eye(s,l);  % Random matrix in [-0.5, 0.5]
       case "sinrand"
            DELTA = rand()*sin(k*pi/4) * eye(s,l);  % Random matrix in [-0.5, 0.5]
       case "aaa"
            DELTA = [sin(k*pi/4) 0
            0 sin(k*pi/4)] ; % Random matrix in [-0.5, 0.5]
         otherwise
            error('Unknown option for CalcDelta: %s', opts);
    end
end