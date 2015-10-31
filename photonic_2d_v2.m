import com.comsol.model.*
import com.comsol.model.util.*

SQ_SIDE_LEN = 375; % in nm
CIRCLE_RADIUS = 0.45; % of SQ_SIDE_LEN
NO_EIGENS = 5;
K_STEP_SIZE = 0.1;

model = ModelUtil.create('Model');

model.modelPath('/your/own/path');

model.modelNode.create('comp1');

model.param.set('a1x', strcat(num2str(SQ_SIDE_LEN, 15), '[nm]')); % parameters
model.param.set('a1y', '0[nm]');
model.param.set('a2x', '0[nm]');
model.param.set('a2y', strcat(num2str(SQ_SIDE_LEN, 15), '[nm]'));
model.param.set('k', 'pi/a1x');
model.param.set('k1', '1');
model.param.set('k2', '1');

model.variable.create('var1'); % variables
model.variable('var1').set('kx', 'k*k1');
model.variable('var1').set('ky', 'k*k2');

model.geom.create('geom1', 2); % geometry
model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature('r1').setIndex('size', 'a1x', 0);
model.geom('geom1').feature('r1').setIndex('size', 'a2y', 1);
model.geom('geom1').feature('r1').set('base', 'center');
model.geom('geom1').run('r1');
model.geom('geom1').feature.create('c1', 'Circle');
model.geom('geom1').feature('c1').set('r', strcat(num2str(CIRCLE_RADIUS, 15), '*a1x'));
model.geom('geom1').run;

model.physics.create('emw', 'ElectromagneticWaves', 'geom1'); % physics
model.physics('emw').feature.create('pc1', 'PeriodicCondition', 1);
model.physics('emw').feature('pc1').selection.set([1 4]);
model.physics('emw').feature('pc1').set('PeriodicType', 1, 'Floquet');
model.physics('emw').feature('pc1').set('kFloquet', {'kx' 'ky' '0'});
model.physics('emw').feature.create('pc2', 'PeriodicCondition', 1);
model.physics('emw').feature('pc2').selection.set([2 3]);
model.physics('emw').feature('pc2').set('PeriodicType', 1, 'Floquet');
model.physics('emw').feature('pc2').set('kFloquet', {'kx' 'ky' '0'});
% the following line wasted my whole Sat. Gotcha!
model.physics('emw').prop('components').set('components', 1, 'outofplane');

model.material.create('mat1'); % material
model.material('mat1').selection.set([2]);
model.material('mat1').name('Air');
model.material('mat1').set('family', 'air');
model.material('mat1').propertyGroup('def').set('relpermeability', '1');
model.material('mat1').propertyGroup('def').set('relpermittivity', '1');
model.material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
model.material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.material('mat1').propertyGroup('def').set('electricconductivity', '0[S/m]');
model.material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
model.material('mat1').propertyGroup('def').set('density', 'rho(pA[1/Pa],T[1/K])[kg/m^3]');
model.material('mat1').propertyGroup('def').set('thermalconductivity', 'k(T[1/K])[W/(m*K)]');
model.material('mat1').propertyGroup('def').set('soundspeed', 'cs(T[1/K])[m/s]');
model.material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.material('mat1').propertyGroup('def').func('eta').set('funcname', 'eta');
model.material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('eta').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.material('mat1').propertyGroup('def').func('Cp').set('funcname', 'Cp');
model.material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('Cp').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.material('mat1').propertyGroup('def').func('rho').set('funcname', 'rho');
model.material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/8.314/T');
model.material('mat1').propertyGroup('def').func('rho').set('dermethod', 'manual');
model.material('mat1').propertyGroup('def').func('rho').set('argders', {'pA' 'd(pA*0.02897/8.314/T,pA)'; 'T' 'd(pA*0.02897/8.314/T,T)'});
model.material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.material('mat1').propertyGroup('def').func('k').set('funcname', 'k');
model.material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('k').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.material('mat1').propertyGroup('def').func('cs').set('funcname', 'cs');
model.material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*287*T)');
model.material('mat1').propertyGroup('def').func('cs').set('dermethod', 'manual');
model.material('mat1').propertyGroup('def').func('cs').set('argders', {'T' 'd(sqrt(1.4*287*T),T)'});
model.material('mat1').propertyGroup('def').addInput('temperature');
model.material('mat1').propertyGroup('def').addInput('pressure');
model.material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat1').propertyGroup('RefractiveIndex').set('n', '1');
model.material.create('mat2');
model.material('mat2').selection.set([1]);
model.material('mat2').propertyGroup('def').set('relpermittivity', {'13'});
model.material('mat2').propertyGroup('def').set('relpermeability', {'1'});
model.material('mat2').propertyGroup('def').set('electricconductivity', {'0'});

model.mesh.create('mesh1', 'geom1'); % mesh
model.mesh('mesh1').automatic(true);
model.mesh('mesh1').autoMeshSize(3); % 9 extremely coarse - 1 extremely fine

model.study.create('std1'); % study
model.study('std1').feature.create('eig', 'Eigenfrequency');
model.study('std1').feature('eig').activate('emw', true);
model.study('std1').feature('eig').set('neigs', num2str(NO_EIGENS));
model.study('std1').feature('eig').set('shift', '1.7e14');
model.sol.create('sol1'); % solution
model.sol('sol1').study('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'eig');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'eig');
model.sol('sol1').feature.create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('shift', '1.7e14');
model.sol('sol1').feature('e1').set('neigs', NO_EIGENS);
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('control', 'eig');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').attach('std1');

%%%%%%%%%%%%%%%%%%%%%%%%%% loop %%%%%%%%%%%%%%%%%%%%%%%%%%

freqs = ones(NO_EIGENS, (1/K_STEP_SIZE+1)*3); % one row one curve
col_idx = 1;

%%% GAMMA - X
% k1 = [0, 1] and k2 = 0
model.param.set('k1', num2str(0, 15));
model.param.set('k2', num2str(0, 15));

k_var = 0;
while k_var <= 1
    model.sol('sol1').runAll; % solve
    sol_info = mphsolinfo(model);
    sol_vals = sol_info.solvals;
    % sanity check
    sol_size = size(sol_vals);
    if sol_size(1) > NO_EIGENS
        %disp('Extra free solutions truncated.');
        sol_vals = sol_vals(1:NO_EIGENS);
    elseif sol_size(1) < NO_EIGENS
        %disp('Not enough solutions returnted. 1 padded.');
        tail = ones(NO_EIGENS-sol_size(1), 1);
        sol_vals = [sol_vals; tail];
    end
    freqs(:, col_idx) = -imag(sol_vals)*SQ_SIDE_LEN*1e-9/(2*pi*3e8);
    
    col_idx = col_idx+1;
    k_var = k_var+K_STEP_SIZE;
    model.param.set('k1', num2str(k_var, 15));
end

%%% X - M
% k1 = 1 and k2 = [0, 1]
model.param.set('k1', num2str(1, 15));
model.param.set('k2', num2str(0, 15));

k_var = 0;
while k_var <= 1
    model.sol('sol1').runAll; % solve
    sol_info = mphsolinfo(model);
    sol_vals = sol_info.solvals;
    % sanity check
    sol_size = size(sol_vals);
    if sol_size(1) > NO_EIGENS
        %disp('Extra free solutions truncated.');
        sol_vals = sol_vals(1:NO_EIGENS);
    elseif sol_size(1) < NO_EIGENS
        %disp('Not enough solutions returnted. 1 padded.');
        tail = ones(NO_EIGENS-sol_size(1), 1);
        sol_vals = [sol_vals; tail];
    end
    freqs(:, col_idx) = -imag(sol_vals)*SQ_SIDE_LEN*1e-9/(2*pi*3e8);
    
    col_idx = col_idx+1;
    k_var = k_var+K_STEP_SIZE;
    model.param.set('k2', num2str(k_var, 15));
end

%%% M - GAMMA
% k1 = [1, 0] and k2 = [1, 0]
model.param.set('k1', num2str(1, 15));
model.param.set('k2', num2str(1, 15));

k_var = 1;
while k_var >= 0
    model.sol('sol1').runAll; % solve
    sol_info = mphsolinfo(model);
    sol_vals = sol_info.solvals;
    % sanity check
    sol_size = size(sol_vals);
    if sol_size(1) > NO_EIGENS
        %disp('Extra free solutions truncated.');
        sol_vals = sol_vals(1:NO_EIGENS);
    elseif sol_size(1) < NO_EIGENS
        %disp('Not enough solutions returnted. 1 padded.');
        tail = ones(NO_EIGENS-sol_size(1), 1);
        sol_vals = [sol_vals; tail];
    end
    freqs(:, col_idx) = -imag(sol_vals)*SQ_SIDE_LEN*1e-9/(2*pi*3e8);
    
    col_idx = col_idx+1;
    k_var = k_var-K_STEP_SIZE;
    model.param.set('k1', num2str(k_var, 15));
    model.param.set('k2', num2str(k_var, 15));
end

plot(0:1:(1/K_STEP_SIZE+1)*3-1, freqs(1, :), 'blue', 'LineWidth', 1.5);
hold on;
plot(0:1:(1/K_STEP_SIZE+1)*3-1, freqs(2, :), 'black', 'LineWidth', 1.5);
plot(0:1:(1/K_STEP_SIZE+1)*3-1, freqs(3, :), 'red', 'LineWidth', 1.5);
plot(0:1:(1/K_STEP_SIZE+1)*3-1, freqs(4, :), 'green', 'LineWidth', 1.5);
plot(0:1:(1/K_STEP_SIZE+1)*3-1, freqs(5, :), 'yellow', 'LineWidth', 1.5);
hold off;
