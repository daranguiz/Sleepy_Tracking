function saver( name )
%SAVER Saves the needed info for the T matrix and it's generation.

if(nargin == 0)
    save('100state25sensRand','T','statesx','statesy','sensx','sensy','AWAKE_COST');
else
    save(name,'T','statesx','statesy','sensx','sensy','AWAKE_COST');
end

end