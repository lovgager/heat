classdef heatEquationWithExchange < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        HeatConductionUIFigure  matlab.ui.Figure
        UIAxes                  matlab.ui.control.UIAxes
        tLabel                  matlab.ui.control.Label
        tSlider                 matlab.ui.control.Slider
        tEditField              matlab.ui.control.NumericEditField
        TEditFieldLabel         matlab.ui.control.Label
        TEditField              matlab.ui.control.NumericEditField
        Panel_2                 matlab.ui.container.Panel
        NEditFieldLabel         matlab.ui.control.Label
        NEditField              matlab.ui.control.NumericEditField
        MEditFieldLabel         matlab.ui.control.Label
        MEditField              matlab.ui.control.NumericEditField
        ButtonGroupPhiPsi       matlab.ui.container.ButtonGroup
        psiEditField            matlab.ui.control.EditField
        Label                   matlab.ui.control.Label
        aEditFieldLabel         matlab.ui.control.Label
        aEditField              matlab.ui.control.NumericEditField
        LEditFieldLabel         matlab.ui.control.Label
        LEditField              matlab.ui.control.NumericEditField
        alphaSlider             matlab.ui.control.Slider
        Label_2                 matlab.ui.control.Label
        alphaEditField          matlab.ui.control.EditField
        cEditField              matlab.ui.control.EditField
        cLabel                  matlab.ui.control.Label
        BuildButton             matlab.ui.control.Button
        NeumannsumSliderLabel   matlab.ui.control.Label
        NeumannsumSlider        matlab.ui.control.Slider
    end

    
    properties (Access = private)
        phi
        phiTest
        psiTest
        psi
        psiPlot
        numericPlot
        c
        cPlot
        t = 0
        numericX
        numericY
        Nsteps = 100
        Msteps = 100
        T = 2
        neumann = 0
        alpha = 0.5
    end
    
    methods (Access = private)
        
        function buildNumeric(app)
            h = pi/app.Nsteps;
            xLocal = 0:h:pi;
            tau = app.t/app.Msteps;
            sigma = tau/h^2;
            y = app.phi;
            
            a = -sigma; 
            b = zeros(1, app.Nsteps-1);
            for i = 1:app.Nsteps
                b(i) = 1 + 2*sigma + tau*app.c(i*h);
            end
            
            for k = 1:app.Msteps
                y = tridiagonalCoef(a, b, y);
            end
            app.numericX = xLocal;
            app.numericY = y;
            plot(app.UIAxes, app.numericX, y, '.', 'LineWidth', 2, 'Color', 'red');     
        end
        
        function buildNumericPhi(app)
            h = pi/app.Nsteps;
            tau = app.T/app.Msteps;
            sigma = tau/h^2;
            y = zeros(1, app.Nsteps - 1);
            for i = 1:app.Nsteps - 1
                y(i) = app.psi(i*h);
            end
            
            a = -sigma; 
            b = zeros(1, app.Nsteps-1);
            for i = 1:app.Nsteps-1
                b(i) = 1 + 2*sigma + tau*app.c(i*h);
            end
            
            app.phi = [0 y 0];
            for periods = 1:app.neumann
                for k = 1:app.Msteps
                    y = tridiagonalCoef(a, b, y);
                end
                app.phi = app.phi + (1 - 1 / app.alpha)^periods * [0 y 0];
            end
            app.phi = app.phi / app.alpha;
        end
        
        function fs = convertStr(~, s)
            fs = '';
            k = 1;
            for i = 1:length(s)
                if ((s(i) == '*') || (s(i) == '^') || (s(i) == '/')) 
                    fs(k) = '.';
                    fs(k + 1) = s(i);
                    k = k + 2;
                else
                    fs(k) = s(i);
                    k = k + 1;
                end
            end
        end
        
        function refreshPlots(app)
            plotPsiC(app);
            buildNumericPhi(app);
            %buildPartSum(app);
            buildNumeric(app);
            legend(app.UIAxes, 'ÿ(x)', 'c(x)', 'Iterative', 'Location', 'best');
        end
        
        function plotPsiC(app)
            hold(app.UIAxes, 'off');
            app.psiPlot = fplot(app.UIAxes, app.psi, [0 pi], 'LineWidth', 2, 'Color', [1 0.56 0.12]);
            hold(app.UIAxes, 'on');
            app.cPlot = fplot(app.UIAxes, app.c, [0 pi], 'LineWidth', 2, 'Color', 'green');
            legend(app.UIAxes, 'ÿ(x)', 'c(x)', 'Location', 'best');
        end
        
        function buildPartSum(app)
            %c(x) = 3/(4x^2), bessel function of first order
            %psi(x) = sqrt(x)*(besselj(1, mu_1 * x/pi) + besselj(1, mu_2 * x/pi)
            bzeros = besselzero(1, 2, 1);
            mu_1 = bzeros(1);
            mu_2 = bzeros(2);
            
            x = linspace(0, pi, 1000);
            res = zeros(1000);
            X_1 = sqrt(x).*besselj(1, x*mu_1/pi);
            X_2 = sqrt(x).*besselj(1, x*mu_2/pi);
            lambda_1 = (mu_1/pi)^2;
            lambda_2 = (mu_2/pi)^2;
            A_1 = 1/(app.alpha + (1-app.alpha)*exp(-lambda_1*app.T));
            A_2 = 1/(app.alpha + (1-app.alpha)*exp(-lambda_2*app.T));
            
            u = A_1*exp(-lambda_1*app.t)*X_1 + A_2*exp(-lambda_2*app.t)*X_2;
            plot(app.UIAxes, x, u, 'LineWidth', 2, 'Color', [0 0.45 0.74]);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            xlim(app.UIAxes, [0 pi]);
            ylim(app.UIAxes, [-inf inf]);
            grid(app.UIAxes, 'on');
        end

        % Value changed function: psiEditField
        function psiEditFieldValueChanged(app, event)
            psiString = convertStr(app, app.psiEditField.Value);
            app.psi = str2func(['@(x)' psiString]);
        end

        % Value changing function: tSlider
        function tSliderValueChanging(app, event)
            app.t = round(event.Value, 3);
            app.tEditField.Value = app.t;
            refreshPlots(app);
        end

        % Value changed function: TEditField
        function TEditFieldValueChanged(app, event)
            app.T = app.TEditField.Value;
            app.tSlider.Limits = [0 app.T];
        end

        % Value changed function: NEditField
        function NEditFieldValueChanged(app, event)
            app.Nsteps = app.NEditField.Value;
        end

        % Value changed function: MEditField
        function MEditFieldValueChanged(app, event)
            app.Msteps = app.MEditField.Value;
            refreshPlots(app);
        end

        % Value changed function: NeumannsumSlider
        function NeumannsumSliderValueChanged(app, event)
            app.neumann = round(app.NeumannsumSlider.Value);
            refreshPlots(app);
        end

        % Value changing function: alphaSlider
        function alphaSliderValueChanging(app, event)
            app.alpha = event.Value;
            app.alphaEditField.Value = sprintf("%.4f", app.alpha);
            refreshPlots(app);
        end

        % Value changed function: cEditField
        function cEditFieldValueChanged(app, event)
            cString = convertStr(app, app.cEditField.Value);
            app.c = str2func(['@(x)' cString]);
        end

        % Button pushed function: BuildButton
        function BuildButtonPushed(app, event)
            app.phiTest = str2func(['@(x)' '(pi/2)^10-(x-pi/2)^10']);
            plotPsiC(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create HeatConductionUIFigure and hide until all components are created
            app.HeatConductionUIFigure = uifigure('Visible', 'off');
            app.HeatConductionUIFigure.Position = [100 100 716 518];
            app.HeatConductionUIFigure.Name = 'Heat Conduction';

            % Create UIAxes
            app.UIAxes = uiaxes(app.HeatConductionUIFigure);
            title(app.UIAxes, '')
            xlabel(app.UIAxes, '')
            ylabel(app.UIAxes, '')
            app.UIAxes.Position = [167 101 499 410];

            % Create tLabel
            app.tLabel = uilabel(app.HeatConductionUIFigure);
            app.tLabel.HorizontalAlignment = 'right';
            app.tLabel.Position = [110 79 25 22];
            app.tLabel.Text = 't =';

            % Create tSlider
            app.tSlider = uislider(app.HeatConductionUIFigure);
            app.tSlider.Limits = [0 2];
            app.tSlider.ValueChangingFcn = createCallbackFcn(app, @tSliderValueChanging, true);
            app.tSlider.Position = [193 93 477 3];

            % Create tEditField
            app.tEditField = uieditfield(app.HeatConductionUIFigure, 'numeric');
            app.tEditField.Editable = 'off';
            app.tEditField.Position = [140 78 38 22];

            % Create TEditFieldLabel
            app.TEditFieldLabel = uilabel(app.HeatConductionUIFigure);
            app.TEditFieldLabel.HorizontalAlignment = 'right';
            app.TEditFieldLabel.Position = [37 78 25 22];
            app.TEditFieldLabel.Text = 'T =';

            % Create TEditField
            app.TEditField = uieditfield(app.HeatConductionUIFigure, 'numeric');
            app.TEditField.Limits = [0 30];
            app.TEditField.ValueChangedFcn = createCallbackFcn(app, @TEditFieldValueChanged, true);
            app.TEditField.Position = [66 78 37 22];
            app.TEditField.Value = 2;

            % Create Panel_2
            app.Panel_2 = uipanel(app.HeatConductionUIFigure);
            app.Panel_2.Position = [9 212 149 42];

            % Create NEditFieldLabel
            app.NEditFieldLabel = uilabel(app.Panel_2);
            app.NEditFieldLabel.HorizontalAlignment = 'right';
            app.NEditFieldLabel.Position = [14 10 19 22];
            app.NEditFieldLabel.Text = 'N =';

            % Create NEditField
            app.NEditField = uieditfield(app.Panel_2, 'numeric');
            app.NEditField.ValueChangedFcn = createCallbackFcn(app, @NEditFieldValueChanged, true);
            app.NEditField.Position = [39 10 31 22];
            app.NEditField.Value = 100;

            % Create MEditFieldLabel
            app.MEditFieldLabel = uilabel(app.Panel_2);
            app.MEditFieldLabel.HorizontalAlignment = 'right';
            app.MEditFieldLabel.Position = [75 10 26 22];
            app.MEditFieldLabel.Text = 'M =';

            % Create MEditField
            app.MEditField = uieditfield(app.Panel_2, 'numeric');
            app.MEditField.ValueChangedFcn = createCallbackFcn(app, @MEditFieldValueChanged, true);
            app.MEditField.Position = [106 10 31 22];
            app.MEditField.Value = 100;

            % Create ButtonGroupPhiPsi
            app.ButtonGroupPhiPsi = uibuttongroup(app.HeatConductionUIFigure);
            app.ButtonGroupPhiPsi.Position = [9 261 149 229];

            % Create psiEditField
            app.psiEditField = uieditfield(app.ButtonGroupPhiPsi, 'text');
            app.psiEditField.ValueChangedFcn = createCallbackFcn(app, @psiEditFieldValueChanged, true);
            app.psiEditField.Position = [35 128 102 22];

            % Create Label
            app.Label = uilabel(app.ButtonGroupPhiPsi);
            app.Label.FontSize = 14;
            app.Label.Position = [9 128 28 22];
            app.Label.Text = 'ÿ =';

            % Create aEditFieldLabel
            app.aEditFieldLabel = uilabel(app.ButtonGroupPhiPsi);
            app.aEditFieldLabel.HorizontalAlignment = 'right';
            app.aEditFieldLabel.FontSize = 14;
            app.aEditFieldLabel.Position = [4 193 25 22];
            app.aEditFieldLabel.Text = 'a =';

            % Create aEditField
            app.aEditField = uieditfield(app.ButtonGroupPhiPsi, 'numeric');
            app.aEditField.Editable = 'off';
            app.aEditField.Position = [35 193 26 22];
            app.aEditField.Value = 1;

            % Create LEditFieldLabel
            app.LEditFieldLabel = uilabel(app.ButtonGroupPhiPsi);
            app.LEditFieldLabel.HorizontalAlignment = 'right';
            app.LEditFieldLabel.FontSize = 14;
            app.LEditFieldLabel.Position = [68 193 25 22];
            app.LEditFieldLabel.Text = 'L =';

            % Create LEditField
            app.LEditField = uieditfield(app.ButtonGroupPhiPsi, 'numeric');
            app.LEditField.Editable = 'off';
            app.LEditField.Position = [100 193 36 22];
            app.LEditField.Value = 3.14159265358979;

            % Create alphaSlider
            app.alphaSlider = uislider(app.ButtonGroupPhiPsi);
            app.alphaSlider.Limits = [0 1];
            app.alphaSlider.ValueChangingFcn = createCallbackFcn(app, @alphaSliderValueChanging, true);
            app.alphaSlider.Position = [12 81 122 3];
            app.alphaSlider.Value = 0.5;

            % Create Label_2
            app.Label_2 = uilabel(app.ButtonGroupPhiPsi);
            app.Label_2.HorizontalAlignment = 'right';
            app.Label_2.FontSize = 14;
            app.Label_2.Position = [39 98 30 22];
            app.Label_2.Text = 'ÿ = ';

            % Create alphaEditField
            app.alphaEditField = uieditfield(app.ButtonGroupPhiPsi, 'text');
            app.alphaEditField.Editable = 'off';
            app.alphaEditField.Position = [74 98 40 22];
            app.alphaEditField.Value = '0.5';

            % Create cEditField
            app.cEditField = uieditfield(app.ButtonGroupPhiPsi, 'text');
            app.cEditField.ValueChangedFcn = createCallbackFcn(app, @cEditFieldValueChanged, true);
            app.cEditField.Position = [36 160 101 22];

            % Create cLabel
            app.cLabel = uilabel(app.ButtonGroupPhiPsi);
            app.cLabel.FontSize = 14;
            app.cLabel.Position = [10 160 25 22];
            app.cLabel.Text = 'c =';

            % Create BuildButton
            app.BuildButton = uibutton(app.ButtonGroupPhiPsi, 'push');
            app.BuildButton.ButtonPushedFcn = createCallbackFcn(app, @BuildButtonPushed, true);
            app.BuildButton.Position = [23 14 100 22];
            app.BuildButton.Text = 'Build';

            % Create NeumannsumSliderLabel
            app.NeumannsumSliderLabel = uilabel(app.HeatConductionUIFigure);
            app.NeumannsumSliderLabel.HorizontalAlignment = 'right';
            app.NeumannsumSliderLabel.Position = [96 42 84 22];
            app.NeumannsumSliderLabel.Text = 'Neumann sum';

            % Create NeumannsumSlider
            app.NeumannsumSlider = uislider(app.HeatConductionUIFigure);
            app.NeumannsumSlider.Limits = [0 10];
            app.NeumannsumSlider.MajorTicks = [0 1 2 3 4 5 6 7 8 9 10];
            app.NeumannsumSlider.MajorTickLabels = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
            app.NeumannsumSlider.ValueChangedFcn = createCallbackFcn(app, @NeumannsumSliderValueChanged, true);
            app.NeumannsumSlider.Position = [192 52 478 3];

            % Show the figure after all components are created
            app.HeatConductionUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = heatEquationWithExchange

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.HeatConductionUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.HeatConductionUIFigure)
        end
    end
end