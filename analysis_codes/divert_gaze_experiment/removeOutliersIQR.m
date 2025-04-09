function cleanedArray = removeOutliersIQR(dataArray)
    % Calculate Q1 and Q3
    Q1 = quantile(dataArray, 0.25);
    Q3 = quantile(dataArray, 0.75);
    
    % Calculate Interquartile Range (IQR)
    IQR = Q3 - Q1;
    
    % Define lower and upper bounds for outliers
    lowerBound = Q1 - 1.5 * IQR;
    upperBound = Q3 + 1.5 * IQR;
    
    % Find outliers
    outliers = (dataArray < lowerBound) | (dataArray > upperBound);
    
    % Remove outliers
    cleanedArray = dataArray(~outliers);
end