%Calculate the normal to the camber line

function norm = calc_norm(x_cam, y_cam)

s = zeros(1, length(x_cam));
grad = zeros(length(x_cam), 2);
norm = zeros(length(x_cam), 2);
for i=2:length(x_cam);
    s(i) = s(i-1)+((x_cam(i)-x_cam(i-1))^2+(y_cam(i)-y_cam(i-1))^2)^0.5;
end
grad(1:end-1, 1) = (x_cam(2:end)-x_cam(1:end-1))./(s(2:end)-s(1:end-1));
grad(1:end-1, 2) = (y_cam(2:end)-y_cam(1:end-1))./(s(2:end)-s(1:end-1));
grad(end, :) = grad(end-1, :);
norm(:, 1) = -grad(:, 2);
norm(:, 2) = grad(:, 1);

end

