// GDAL.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "gdal_priv.h"
#pragma comment( lib, "opencv_highgui249d.lib")
#pragma comment( lib, "opencv_core249d.lib")
 
#include "ASiftDetector.h"
#include "utils.h"

// ��ȡsift����
Mat ExtractSIFTFeature(const string &imgfn, vector<KeyPoint> &points)  {
    auto img = imread(imgfn, true);    //imgfn: image file name
    SiftFeatureDetector detector;
    vector<KeyPoint> keypoints;
    detector.detect(img, keypoints);
    SiftDescriptorExtractor extractor;
    Mat descriptors;
    if (!keypoints.size()) {
        return Mat();
    }
    // extract sift feature
    extractor.compute(img, keypoints, descriptors);
    //points.resize(keypoints.size());
    for (int i = 0; i < keypoints.size(); i++){
        points.push_back(keypoints[i]);
    }
    return descriptors;
}

/*
    //ע���ļ���ʽ
    GDALAllRegister();
   
    const char* pszFile = "C:\\Users\\Xyy\\Desktop\\SIFT\\GF_L2A.tif";
    GDALDataset *poDataset;
    //ʹ��ֻ����ʽ��ͼ��
    poDataset = (GDALDataset*) GDALOpen( pszFile,GA_ReadOnly );
    if( poDataset == NULL )
    {
        printf( "File: %s���ܴ򿪣�\n",pszFile);
        return 0;
    }
 
    //���ͼ��ĸ�ʽ��Ϣ
    printf( "Driver:%s/%s\n",
        poDataset->GetDriver()->GetDescription(),
        poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME) );
 
    //���ͼ��Ĵ�С�Ͳ��θ���
    printf( "Size is%dx%dx%d\n",
        poDataset->GetRasterXSize(),poDataset->GetRasterYSize(),
        poDataset->GetRasterCount());
 
    //���ͼ���ͶӰ��Ϣ
    if( poDataset->GetProjectionRef() != NULL )
        printf( "Projectionis `%s'\n", poDataset->GetProjectionRef() );
 
    //���ͼ�������ͷֱ�����Ϣ
    double adfGeoTransform[6];
    if( poDataset->GetGeoTransform( adfGeoTransform) == CE_None )
    {
        printf( "Origin =(%.6f,%.6f)\n",
            adfGeoTransform[0], adfGeoTransform[3]);
 
        printf( "PixelSize = (%.6f,%.6f)\n",
            adfGeoTransform[1], adfGeoTransform[5]);
    }
 
    GDALRasterBand *poBand;
    int            nBlockXSize, nBlockYSize;
    int            bGotMin, bGotMax;
    double         adfMinMax[2];
 
    //��ȡ��һ������
    poBand = poDataset->GetRasterBand( 1 );
 
    //��ȡͼ��Ŀ��С�����
    poBand->GetBlockSize(&nBlockXSize, &nBlockYSize );
    printf( "Block=%dx%dType=%s, ColorInterp=%s\n",
        nBlockXSize, nBlockYSize,
        GDALGetDataTypeName(poBand->GetRasterDataType()),
        GDALGetColorInterpretationName(
        poBand->GetColorInterpretation()));
 
    //��ȡ�ò��ε����ֵ��Сֵ�������ȡʧ�ܣ������ͳ��
    adfMinMax[0] = poBand->GetMinimum( &bGotMin);
    adfMinMax[1] = poBand->GetMaximum( &bGotMax);
 
    if( ! (bGotMin&& bGotMax) )
        GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
 
    printf( "Min=%.3fd,Max=%.3f\n", adfMinMax[0], adfMinMax[1] );
 
    //���ͼ��Ľ�������Ϣ
    if( poBand->GetOverviewCount() > 0 )
        printf( "Band has%d overviews.\n", poBand->GetOverviewCount() );
 
    //���ͼ�����ɫ����Ϣ
    if( poBand->GetColorTable() != NULL)
        printf( "Band hasa color table with %d entries.\n",
        poBand->GetColorTable()->GetColorEntryCount() );
 
    float *pafScanline;
    int   nXSize = poBand->GetXSize();
   
    //��ȡͼ��ĵ�һ������
    pafScanline = (float*) CPLMalloc(sizeof(float)*nXSize);
    poBand->RasterIO(GF_Read, 0, 0, nXSize,1, 
        pafScanline, nXSize,1, GDT_Float32, 0, 0 );
 
    CPLFree(pafScanline);
 
    //�ر��ļ�
    GDALClose((GDALDatasetH)poDataset);
	
	*/



int main(int argc, const char * argv[]) {
    // insert code here...
    string imgfn = "C:\\Users\\Xyy\\Desktop\\SIFT\\GF01_MS2_cut.tif";
    Mat queryImage, queryBlackImg, qDescriptor;
    queryImage = imread(imgfn);
    vector<KeyPoint> qKeypoints;
    qDescriptor = ExtractSIFTFeature(imgfn, qKeypoints);
    /*
	drawKeypoints(queryImage, qKeypoints, queryBlackImg);
    imshow("Sift", queryBlackImg);
    cvWaitKey(0);
	*/
    
    string objFileName = "C:\\Users\\Xyy\\Desktop\\SIFT\\15.tif";
    Mat objectImage, objBlackImg, objDesriptor;
    objectImage = imread(objFileName);
    vector<KeyPoint> objKeypoints;
    objDesriptor = ExtractSIFTFeature(objFileName, objKeypoints);
    /*
	drawKeypoints(objectImage, objKeypoints, objBlackImg);
    imshow("Sift", objBlackImg);
    cvWaitKey(0);
	*/
    
    //ֱ����1���ڣ�SIFT��ƥ���
    FlannBasedMatcher matcher;
	
    vector< DMatch > matches1;
    matcher.match(qDescriptor, objDesriptor, matches1);
    findInliers(qKeypoints, objKeypoints, matches1, imgfn, objFileName); // ͨ��1����ƥ��ĵ��
    
    // ʹ��Lowe��raw feature����ƥ��
    vector<Point2f> qeK, obK;
    vector<vector<DMatch>> matches;
    vector<DMatch> good_matches2;
    matcher.knnMatch(qDescriptor, objDesriptor, matches, 2);
    for (size_t i = 0; i < matches.size(); i++){
        if (matches[i][0].distance < 0.8*matches[i][1].distance){
            good_matches2.push_back(matches[i][0]);
            qeK.push_back(qKeypoints[matches[i][0].queryIdx].pt);
            obK.push_back(objKeypoints[matches[i][1].trainIdx].pt);
        }
    }
    findInliers(qKeypoints, objKeypoints, good_matches2, imgfn, objFileName); // ͨ��Loweƥ��ĵ��

    
    /*
	ASiftDetector asiftDetector;
    vector<KeyPoint> asiftKeypoints_query;
    Mat asiftDescriptors_query;
    asiftDetector.detectAndCompute(queryImage, asiftKeypoints_query, asiftDescriptors_query);
    vector<KeyPoint> asiftKeypoints_object;
    Mat asiftDescriptors_object;
    asiftDetector.detectAndCompute(objectImage, asiftKeypoints_object, asiftDescriptors_object);*/
    
    /*Matching descriptor vectors using FLANN matcher, ASIFT��ƥ���
    std::vector< DMatch > asiftMatches;
    matcher.match(asiftDescriptors_query, asiftDescriptors_object, asiftMatches);
    findInliers(asiftKeypoints_query, asiftKeypoints_object, asiftMatches, imgfn, objFileName);*/
    
    // ʹ�����ú�����ƥ����  Mat img_matches;
     /*
	 drawMatches( queryImage, qKeypoints, objectImage, objKeypoints,  matches, img_matches, Scalar::all(-1), Scalar::all(-1),
     vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
     imshow("Good Matches & Object detection", img_matches);
     waitKey(0);
	 */
    
    return 0;
}

	
	
		
		   
	/*vector<Point2f> qeK, obK;
    vector<vector<DMatch>> matches;
    vector<DMatch> good_matches2;
    matcher.knnMatch(qDescriptor, objDesriptor, matches, 2);
    for (size_t i = 0; i < matches.size(); i++){
        if (matches[i][0].distance < 0.8*matches[i][1].distance){
            good_matches2.push_back(matches[i][0]);
            qeK.push_back(result[matches[i][0].queryIdx].pt);

			//qeK.push_back(qKeypoints[matches[i][0].queryIdx].pt);

            obK.push_back(objKeypoints[matches[i][1].trainIdx].pt);
        }
    }
    findInliers(result, objKeypoints, good_matches2, imgfn, objFileName); // ͨ��Loweƥ��ĵ��
	*/
	/*findInliers(qKeypoints, objKeypoints, good_matches2, imgfn, objFileName); // ͨ��Loweƥ��ĵ��
		vector< DMatch > matches1;
		Mat testdescriptors;
		SiftFeatureDetector detector;
		SiftDescriptorExtractor extractor;
		extractor.compute(queryImage, result, testdescriptors);
		matcher.match(testdescriptors, objDesriptor, matches1);
		*/
    	/*
    vector< DMatch > matches1;
    matcher.match(qDescriptor, objDesriptor, matches1);
    findInliers(qKeypoints, objKeypoints, matches1, imgfn, objFileName); // ͨ��1����ƥ��ĵ��
    */
    /*ASiftDetector asiftDetector;
    vector<KeyPoint> asiftKeypoints_query;
    Mat asiftDescriptors_query;
    asiftDetector.detectAndCompute(queryImage, asiftKeypoints_query, asiftDescriptors_query);
    vector<KeyPoint> asiftKeypoints_object;
    Mat asiftDescriptors_object;
    asiftDetector.detectAndCompute(objectImage, asiftKeypoints_object, asiftDescriptors_object);*/
    
    /*Matching descriptor vectors using FLANN matcher, ASIFT��ƥ���
    std::vector< DMatch > asiftMatches;
    matcher.match(asiftDescriptors_query, asiftDescriptors_object, asiftMatches);
    findInliers(asiftKeypoints_query, asiftKeypoints_object, asiftMatches, imgfn, objFileName);*/
    
    // ʹ�����ú�����ƥ����  Mat img_matches;
     /*
	 drawMatches( queryImage, qKeypoints, objectImage, objKeypoints,  matches, img_matches, Scalar::all(-1), Scalar::all(-1),
     vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS );
     imshow("Good Matches & Object detection", img_matches);
     waitKey(0);
	 */
	 /*drawKeypoints(queryImage, qKeypoints, queryBlackImg);
    imshow("Sift", queryBlackImg);
    cvWaitKey(0);*/
		