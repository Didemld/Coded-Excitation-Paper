function [lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(calData)

    % Find the maximum point in the data (assuming PSF center)
    [Max_r, Max_c] = find(calData == max(calData(:)));

    % Extract lateral PSF (row-based)
    lateralPSF = squeeze(calData(Max_r, :));
    lateralPSF = mag2db(abs(lateralPSF));
    lateralPSF = lateralPSF - max(lateralPSF);

    % Extract axial PSF (column-based)
    axialPSF = squeeze(calData(:, Max_c));
    axialPSF = mag2db(abs(axialPSF));
    axialPSF = axialPSF - max(axialPSF);

    % Compute FWHM for lateral PSF
    [maxIndex_L] = find(lateralPSF == max(lateralPSF)); % Find maximum PSF
    halfMax_L = lateralPSF(maxIndex_L) - 6; % Calculate half maximum for lateral PSF

    % Find left and right side half
    leftIndex_L = find(lateralPSF(1:maxIndex_L) < halfMax_L, 1, 'last');
    if isempty(leftIndex_L)
        leftIndex_L = 1; % Use terminal point if not found
    end

    rightIndex_L = find(lateralPSF(maxIndex_L:end) < halfMax_L, 1, 'first') + maxIndex_L - 1;
    if isempty(rightIndex_L)
        rightIndex_L = length(lateralPSF); % Use terminal point if not found
    end

    % Calculate FWHM for lateral PSF
    FWHM_L = 0.025*((rightIndex_L) - (leftIndex_L));

    % Compute FWHM for axial PSF
    [maxIndex_A] = find(axialPSF == max(axialPSF)); % Find maximum PSF
    halfMax_A = axialPSF(maxIndex_A) - 6; % Calculate half maximum for axial PSF

    % Find left and right side half
    leftIndex_A = find(axialPSF(1:maxIndex_A) < halfMax_A, 1, 'last');
    if isempty(leftIndex_A)
        leftIndex_A = 1; % Use terminal point if not found
    end

    rightIndex_A = find(axialPSF(maxIndex_A:end) < halfMax_A, 1, 'first') + maxIndex_A - 1;
    if isempty(rightIndex_A)
        rightIndex_A = length(axialPSF); % Use terminal point if not found
    end

    % Calculate FWHM for axial PSF
    FWHM_A = 0.025*((rightIndex_A) - (leftIndex_A));

end