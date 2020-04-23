function [Qbar] = degradeQbar(Qbar, failed_plies)
%DEGRADEQBAR Fully degrades stiffness of failed plies in a laminate
%   Detailed explanation goes here

n_failed_plies = length(failed_plies);
for i=1:n_failed_plies
    idx = failed_plies(i);
    Qbar(:,:, idx) = zeros(3);
end

end

