//
//  main.cpp
//  exx
//
//  Created by Kuanhsun wu on 2018/12/15.
//  Copyright © 2018 Kuanhsun wu. All rights reserved.
//
#include <stdio.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <dirent.h> // for linux systems
#include <sys/stat.h> // for linux systems
#include <opencv2/ml/ml.inl.hpp>
#include "opencv2/core/utility.hpp"
#include <map>
#include <opencv2/dnn.hpp>
#include <opencv2/plot.hpp>
#include <cmath>
#include <fstream>


using namespace cv;
using namespace std;
using namespace cv::ml;
using namespace dnn;


void visualizeHOG(cv::Mat img, std::vector<float> &feats, cv::HOGDescriptor hog_detector, int scale_factor = 3);
void feature_extract(Mat &image,Mat &des);
void readFilenames(std::vector<string> &filenames, const string &directory)
{
#ifdef WINDOWS
    HANDLE dir;
    WIN32_FIND_DATA file_data;
    
    if ((dir = FindFirstFile((directory + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
        return; /* No files found */
    
    do {
        const string file_name = file_data.cFileName;
        const string full_file_name = directory + "/" + file_name;
        const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
        
        if (file_name[0] == '.')
            continue;
        
        if (is_directory)
            continue;
        
        filenames.push_back(full_file_name);
    } while (FindNextFile(dir, &file_data));
    
    FindClose(dir);
#else
    DIR *dir;
    class dirent *ent;
    class stat st;
    
    dir = opendir(directory.c_str());
    while ((ent = readdir(dir)) != NULL) {
        const string file_name = ent->d_name;
        const string full_file_name = directory + "/" + file_name;
        
        if (file_name[0] == '.')
            continue;
        
        if (stat(full_file_name.c_str(), &st) == -1)
            continue;
        
        const bool is_directory = (st.st_mode & S_IFDIR) != 0;
        
        if (is_directory)
            continue;
        
        //        filenames.push_back(full_file_name); // returns full path
        filenames.push_back(file_name); // returns just filename
    }
    closedir(dir);
#endif
} // GetFilesInDirectory
//Mat feats, labels; // start empty
void aug_flip (string path){
    
    for(int j = 0 ; j < 4 ; j++){
        String folder = (path + to_string(j) + "/");
//                cout <<  folder << endl;
        vector<string> filenames;
        vector<float> descriptors;
        readFilenames(filenames, folder);
        //        cout << "number of sample "<< filenames.size() << endl;
        for(int  i = 0; i < filenames.size() ; ++i)
        {
            Mat image = imread(folder + filenames[i]);
            Mat h_flip ;
            
            
            
//            blur(image, h_flip, Size(1,1));
//            cv::imwrite(folder + filenames[i]+"_blur_0.jpg", h_flip);
//            blur(image, h_flip, Size(2,2));
//            cv::imwrite(folder + filenames[i]+"_blur_1.jpg", h_flip);
//            blur(image, h_flip, Size(3,3));
//            cv::imwrite(folder + filenames[i]+"_blur_2.jpg", h_flip);
//            blur(image, h_flip, Size(4,4));
//            cv::imwrite(folder + filenames[i]+"_blur_3.jpg", h_flip);
////
//
//            image.convertTo(h_flip, -1, 0.5, 50);
//            cv::imwrite(folder + filenames[i]+"_illum_0.jpg", h_flip);
//            image.convertTo(h_flip, -1, 1, 50);
//            cv::imwrite(folder + filenames[i]+"_illum_1.jpg", h_flip);
//            image.convertTo(h_flip, -1, 1.5, 50);
//            cv::imwrite(folder + filenames[i]+"_illum_2.jpg", h_flip);
//            image.convertTo(h_flip, -1, 2, 50);
//            cv::imwrite(folder + filenames[i]+"_illum_3.jpg", h_flip);

             resize(image, h_flip, Size(image.size().height * 1.5,image.size().width * 1.5),INTER_CUBIC );
             cv::imwrite(folder + filenames[i]+"_scale_0.jpg", h_flip);
            resize(image, h_flip, Size(image.size().height * 0.5,image.size().width * 0.5),INTER_CUBIC );
            cv::imwrite(folder + filenames[i]+"_scale_1.jpg", h_flip);
            resize(image, h_flip, Size(image.size().height * 2,image.size().width * 2),INTER_CUBIC );
            cv::imwrite(folder + filenames[i]+"_scale_3.jpg", h_flip);


            cv::flip(image, h_flip,0);
            cv::imwrite(folder + filenames[i]+"_flip_0.jpg", h_flip);
            cv::flip(image, h_flip,1);
            cv::imwrite(folder + filenames[i]+"_flip_1.jpg", h_flip);
            cv::flip(image, h_flip,-1);
            cv::imwrite(folder + filenames[i]+"_flip_-1.jpg", h_flip);
//            namedWindow("Display window", WINDOW_AUTOSIZE);
//            imshow("Display window", image);
//            waitKey(0);
//            imshow("Display window", h_flip);
//            waitKey(0);
            /*
            Mat dist;
            resize(image, image, Size(128,128),INTER_CUBIC );
            cv::Point2f center((image.cols-1)/2.0, (image.rows-1)/2.0);
//            cv::Mat rot = cv::getRotationMatrix2D(center, 45, 1.0);
//            warpAffine(image,dist , rot, Size(image.cols, image.rows));
//            cv::imwrite(folder + filenames[i]+"_rotate_0.jpg", dist);
            cv::Mat rot = cv::getRotationMatrix2D(center, 90, 1.0);
            warpAffine(image,dist , rot, Size(image.rows, image.cols));
            cv::imwrite(folder + filenames[i]+"_rotate_1.jpg", dist);
//            rot = cv::getRotationMatrix2D(center, 135, 1.0);
//            warpAffine(image,dist , rot, Size(image.cols, image.rows));
//            cv::imwrite(folder + filenames[i]+"_rotate_2.jpg", dist);
            rot = cv::getRotationMatrix2D(center, 180, 1.0);
            warpAffine(image,dist , rot, Size(image.cols, image.rows));
            cv::imwrite(folder + filenames[i]+"_rotate_3.jpg", dist);
//            rot = cv::getRotationMatrix2D(center, 225, 1.0);
//            warpAffine(image,dist , rot, Size(image.cols, image.rows));
//            cv::imwrite(folder + filenames[i]+"_rotate_4.jpg", dist);
            rot = cv::getRotationMatrix2D(center, 270, 1.0);
            warpAffine(image,dist , rot, Size(image.rows, image.cols));
            cv::imwrite(folder + filenames[i]+"_rotate_5.jpg", dist);
//            rot = cv::getRotationMatrix2D(center, 315, 1.0);
//            warpAffine(image,dist , rot, Size(image.cols, image.rows));
//            cv::imwrite(folder + filenames[i]+"_rotate_6.jpg", dist);
             */
            

        }
    }
}
void data (string path,Mat &feats,Mat &labels){
    
    for(int j = 0 ; j < 4 ; j++){
        String folder = (path + to_string(j) + "/");

        vector<string> filenames;
        readFilenames(filenames, folder);
        
//        cout << "number of sample "<< filenames.size() << endl;
        
        for(int  i = 0; i < filenames.size() ; ++i)
            {
                Mat image = imread(folder + filenames[i]);
                Mat des ;
                feature_extract(image, des);
                feats.push_back(des);
//                cout<<des.size() ;
                
                labels.push_back(j);
//                cout <<  labels.at<Vec3b>(i,0) << endl;
//                cout << labels.size() << endl;
            }
    }
}


void create( cv::Ptr<cv::ml::DTrees> dtree[], int num){
    
    for (int i = 0; i < num; i++){
        
        dtree[i] = cv::ml::DTrees::create();
        dtree[i]->setMaxCategories(6);
        dtree[i]->setMaxDepth(20);//5
        dtree[i]->setMinSampleCount(5);//5
        dtree[i]->setCVFolds(0);
//        dtree[i]->setUseSurrogates(false);
//        dtree[i]->setUse1SERule(true);
//        dtree[i]->setTruncatePrunedTree(true);
        dtree[i]->setPriors(cv::Mat());
        
    }
}
void training ( cv::Ptr<cv::ml::DTrees> dtree[],Mat &train,Mat &train_label,int num){
    for (int i = 0; i < num; i++){
//        cout << train.size().width<<endl ;
//        int residual  = 0 ;
//        if(train.size().width % 2 != 0)
//            residual = 0 ;
//        Mat select_feature = Mat::ones((int)(train.size().width * 0.9 + residual) , 1, CV_8U );
//        Mat unselect_feature = Mat::zeros((int)(train.size().width * 0.1)+1    , 1, CV_8U );
//        select_feature.push_back(unselect_feature);
//        srand(time(0));
//        randShuffle(select_feature);
////
////        cout<< select_feature.size()<< endl ;
////
////
//        cv::Ptr<TrainData> train_data =  ml::TrainData::create(train, ml::ROW_SAMPLE, train_label,select_feature) ;
        cv::Ptr<TrainData> train_data =  ml::TrainData::create(train, ml::ROW_SAMPLE, train_label) ;

//        cout << train_data->getNVars()  << endl ;
//        train_data->setTrainTestSplitRatio(0.7,true);
//        srand(time(0));
//        float n = (1 - 0.5) * rand() / (RAND_MAX + 1.0) + 0.5;

//        cout << n <<endl ;
        train_data->setTrainTestSplitRatio(0.7,true);


 //        cout << train_data->getTrainSampleIdx()  << endl ;

        
        dtree[i]->train(train_data);
        //        printf( "train error: %f\n", dtree[i]->calcError(, false, noArray()) );
//        printf( "test error: %f\n\n", dtree[i]->calcError(train_data, true, noArray()) );
    }
}

void predict (cv::Ptr<cv::ml::DTrees> dtree[], Mat &test,Mat &test_result,Mat &prob, int num){
    Mat result[num], vote = Mat::zeros(4,test.size().height,CV_32F);;

    for (int i = 0; i < num; i++){
        
        dtree[i]->predict(test,result[i]) ;
        
        for(int j = 0 ; j < test.size().height; j ++){
            vote.at<float>(result[i].at<float>(0,j),j)++ ;
//            cout << result[i].at<float>(0,j) <<endl ;


        }
//        cout<<"----------------------------\n";
//        cout<<vote ;
//        cout << result[i] <<endl ;
//        cout << test.size().height << endl ;
    }
  
//        cout << vote.size() << endl ;
    for(int i = 0 ; i < test.size().height ; i++){
        float max = 0 , max_vote = 0 ;
        for(int j = 0 ;  j < 4 ; j++ ){
            
//            cout<<i<<" "<<j << " "<< vote.at<float>(j,i) << endl ;
            if ( vote.at<float>(j,i) > max){
                max = vote.at<float>(j,i) ;
               
                max_vote = j ;
            }
        }
//        cout<<"max "<<max<<" vote "<<max_vote <<endl ;
//        cout<<"-----------------------" <<endl ;
        prob.push_back((float)max/num);
        test_result.push_back(max_vote) ;
//        cout << "max_vote"<<max_vote<< endl ;
    }
    
}
void feature_extract(Mat &image,Mat &des){
   
    vector<float> descriptors;
    
    Mat kernel =  (Mat_<double>(5,5) << -1,-1,-1,-1,-1,
                   -1,2,2,2,-1,
                   -1,2,8,2,-1,
                   -2,2,2,2,-1,
                   -1,-1,-1,-1,-1);
    kernel = kernel / 8.0 ;
//    blur(image, image, Size(10,10));
//    cv::filter2D(image,image , -1, kernel);
//    cv::cvtColor(image, image, cv::COLOR_BGR2GRAY);

//    threshold( image, image, 80,255,THRESH_BINARY );

//                                namedWindow("Display window", WINDOW_AUTOSIZE);
//                                imshow("Displa y window", image);
//                                waitKey(0);
//
    if(!image.data)
        cerr << "Problem loading image!!!" << endl;
    resize(image, image, Size(48,48),INTER_CUBIC );
    //                namedWindow("Display window", WINDOW_AUTOSIZE);
    //                             imshow("Display window", image);
    //                             waitKey(0);
    HOGDescriptor hog(
                      Size(48,48), //winSize
                      Size(24,24), //blocksize
                      Size(12,12), //blockStride,
                      Size(12,12), //cellSize,
                      9, //nbins,
                      1, //derivAper,
                      -1, //winSigma,
                      0, //histogramNormType,
                      0.2, //L2HysThresh,
                      1,//gammal correction,
                      64,//nlevels=64
                      1);//Use signed gradients
    
    hog.compute(image,descriptors,Size(12,12),Size(0,0));

    
    des = Mat(descriptors).clone();
    des.convertTo(des, CV_32F);     // ml needs float data
    
    des = des.reshape(1,1);
//    cout<<des.size()<<endl ;
    
}
struct object {
    Rect box ;
    float  id, prob;
};
bool compare_o(object i1, object i2)
{
    return (i1.prob > i2.prob);
}
int main( int argc, char** argv )
{

    Mat train, train_label, train_result, prob ;
    int num = 60 ;
    string folder = ("./data/task3/train/0");
//    aug_flip(folder);
    
    data( folder,train,train_label);
    
    cout << "Vector of floats via Mat = " << train.size()<<endl ;

   

    cv::Ptr<cv::ml::DTrees> dtree[200];
    create(dtree,num);
    training(dtree,train,train_label,num);
    predict(dtree, train, train_result,prob, num);

    Mat diff_train_t3;
    train_result.convertTo(train_result, CV_32S);

    diff_train_t3 = train_result - train_label ;
//    cv::hconcat(train_result, train_label, train_result) ;
//    cout<<train_result;
    

    cout << "train err = "<< countNonZero(diff_train_t3) << endl;
    
    int TP = 0, rank = 0 ;
    std::vector<double> precision,recal;
   
    for (int l = 0 ;l < 44 ; l ++){
        string S ;
        if (l < 10){
            S= "./data/task3/test/000" +to_string(l) + ".jpg" ;
        }else{
            S = "./data/task3/test/00" +to_string(l) + ".jpg" ;
        }
        
//        cout<< S ;
        Mat A = imread(S);
       
    std::vector<float> confidences, classid;
    std::vector<object> object_d , object_r;
    std::vector<Rect> boxes;
        object tmp_o ;
    int g = 0 ;

    for (int k = 0 ; k < 15; k++){
        int stride_x = 16 , stride_y = 16 ,sample_h = 60 + k * 10  , sample_w = 60 + k * 10;
         for(int i = 0 ; i < A.rows - sample_w + 1  ; i = i  + stride_y){
             for (int j = 0 ; j < A.cols - sample_h + 1  ;j = j + stride_x){
    //             cout << "i = "<< i <<"j = "<< j<<endl ;
                 Mat tmp ;
                 A(Rect(j,i, sample_h , sample_w)).copyTo(tmp);
    //             namedWindow("Display window", WINDOW_AUTOSIZE);
    //             imshow("Display window", tmp);
    //             waitKey(0);
                 Mat des;
                 feature_extract(tmp, des);
    //             cout << des.size() << endl;

                
                 Mat re ,prob_test;
                 predict(dtree, des, re, prob_test, num);
                 if (re.at<float>(0,0) != 3 && prob_test.at<float>(0,0) > 0.5){
                     g++;
//                     cout<<"id "<< re.at<float>(0,0) <<" prob :"<<prob_test.at<float>(0,0)<<endl ;
//
//                                  namedWindow("Display window", WINDOW_AUTOSIZE);
//                                  imshow("Display window", A(Rect(j,i ,  sample_h , sample_w)));
//                                  waitKey(0);
                     boxes.push_back(Rect(j,i ,  sample_h , sample_w));
                     confidences.push_back(prob_test.at<float>(0,0));
                     classid.push_back(re.at<float>(0,0)) ;
                     tmp_o.box = Rect(j,i ,  sample_h , sample_w);
                     tmp_o.prob = prob_test.at<float>(0,0);
                     tmp_o.id = re.at<float>(0,0) ;
                     object_d.push_back(tmp_o);
                     
                    }
                 }
             }
    
    }
//    cout<<"g = "<<g<<endl ;

        
        float thr,iou ;
        int area,xx1,xx2,yy1,yy2, area_h,area_w;
        
        sort(object_d.begin(), object_d.end(), compare_o);

        for(int x = 0 ; x < object_d.size() ; x ++){
            sort(object_d.begin(), object_d.end(), compare_o);
            
            if (object_d[x].prob == 0){
                break;
            }
            object_r.push_back(object_d[x]);
            for(int y = x+1 ;y <  object_d.size() ; y ++){
//                Mat v ;
//                A.copyTo(v) ;
//
//                //        cout << "Intervals sorted by start time : \n";
//                //        for (auto x : object_d)
//                //            cout <<x.prob <<endl;
//
//                rectangle( v, object_d[x].box,  Scalar(255,0,0), 2, 8, 0 );
//
//                rectangle( v, object_d[y].box, Scalar(0,255,0), 2, 8, 0 );
                xx1 = max(object_d[x].box.tl().x, object_d[y].box.tl().x);
                yy1 = max(object_d[x].box.tl().y, object_d[y].box.tl().y);
                xx2 = min(object_d[x].box.br().x, object_d[y].box.br().x);
                yy2 = min(object_d[x].box.br().y, object_d[y].box.br().y);
                area_w = max(0, xx2 - xx1 + 1);
                area_h = max(0, yy2 - yy1 + 1);
                area = area_w * area_h ;
                iou = (float)area / (float)(object_d[x].box.area() + object_d[y].box.area() - area);
//                cout << "iou "<<iou<<endl ;
//                imshow("Display window", v);
//                waitKey(0);

                if (iou > 0.1){
                    object_d[y].prob = 0 ;
                }

                //                cout <<boxes[x].br()<<"  "<<boxes[y].br()<<"  "<<xx1<<"\n";
            }
                
        }
        for (size_t i = 0; i < 3; ++i)
        {

            RNG rng(1235);
            Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );

            rectangle( A, object_r[i].box, color, 2, 8, 0 );
            string title = to_string((int)object_r[i].id) + "  "+ to_string(object_r[i].prob) ;
            putText(A, title, Point( object_r[i].box.x,object_r[i].box.y ), 0, 0.4, Scalar(0,0,255),1);
            
        }
        imshow("Display window", A);
        string z = "my" + to_string(l) + ".jpg" ;
        cv::imwrite(z, A);



        
////////////////////NMSBoxes//////////////////////////////////////////////
            std::vector<int> indices, max_3;
        int max_index =  0 ,max_i = 0;
            NMSBoxes(boxes, confidences, 0, 0.1, indices);
        ///find three bounding boxes with the heightest confidence
        for(int j = 0 ; j < 3 ; j++){
            float max = 0 ;
            for (int i = 0; i < indices.size(); ++i)
            {
        
                int idx = indices[i];
              
                float P = confidences[idx] ;
                if (P > max ){
                    max = P;
                    max_index = idx ;
                    max_i = i;
                }
            }
        max_3.push_back(max_index);
        indices.erase(indices.begin() + max_i);
        }
        
        
        for(int i = 0 ; i < max_3.size() ; i ++){
            int idx = max_3[i] ;
            Rect box = boxes[idx];
            float id = classid[idx] ;
            float P = confidences[idx] ;
            RNG rng(12345);
            Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) );
        
            rectangle( A, box, color, 2, 8, 0 );
            string title = to_string((int)id) + "  "+ to_string(P) ;
            putText(A, title, Point( box.x,box.y ), 0, 0.4, Scalar(1,0,0),1);
          
        }
   
//////////////////////////Evaluate the detection result///////
    string gt ;

        if (l < 10){
            gt= "./data/task3/gt/000" +to_string(l) + ".gt.txt" ;
        }else{
            gt =  "./data/task3/gt/00" +to_string(l) + ".gt.txt" ;
        }
        fstream myfile;
        char line[100];
    
    
        myfile.open (gt,ios::in);
        while(myfile.getline(line,sizeof(line),'\n')){
            char * pch;
            pch = strtok (line," ");
            std::vector<int> ans ;
            while (pch != NULL)
            {
                //cout<< atoi(pch)<<endl ;
                int p = atoi(pch) ;
                pch = strtok (NULL, " ");
                ans.push_back(p);
            }
            rank = rank + 1 ;
            cout<<" rank "<<rank <<endl;
            int count = 0 ;
            for(int i = 0 ; i < 3 ; i ++){
                float iou = 0 ;
                int area,xx1,xx2,yy1,yy2, area_h,area_w;
                
               
                Rect box = object_r[i].box;
                float id = object_r[i].id ;
                if ((int)id == ans[0]){
                    count ++ ;
                    xx1 = max(box.tl().x,ans[1]);
                    yy1 = max(box.tl().y,ans[2]);
                    xx2 = min(box.br().x,ans[3]);
                    yy2 = min(box.br().y,ans[4]);
                    area_w = max(0, xx2 - xx1 + 1);
                    area_h = max(0, yy2 - yy1 + 1);
                    area = area_w * area_h ;
                    iou = (float)area / (float)(box.area() + (ans[3] - ans[1]) * (ans[4] - ans[2]) - area);
                    if(iou > 0.5){
                        TP = TP + 1 ;
                    }
                }
            }
            precision.push_back((float)TP / -rank);
            cout<<" precision "<< (float)TP / rank<<endl;
            recal.push_back((float)TP / (44 * 3));
            cout << " recal "<< (float)TP / (44 * 3)<<endl ;
        }
 
        myfile.close();
     
    imshow("Display window", A);
    //waitKey(-1);
    string w = to_string(l) + ".jpg" ;
    cv::imwrite(w, A);

    }
//arbitrarilly selected three highest bounding box
// so FP = FN
    Mat display ;
    cv::Ptr<cv::plot::Plot2d> plot = cv::plot::Plot2d::create(Mat(recal), Mat(precision) );
    plot->setShowText(true);
    plot->setGridLinesNumber(10);
   
    plot->render(display);
    imshow("Plot", display);
    waitKey(0);
    cv::imwrite("PR curve_thr_0.5.jpg", display);
    
    
    
    
   
    


    
    return 0;
}










void visualizeHOG(cv::Mat img, std::vector<float> &feats, cv::HOGDescriptor hog_detector, int scale_factor ) {

    cv::Mat visual_image;
    resize(img, visual_image, cv::Size(img.cols * scale_factor, img.rows * scale_factor));
//
    int n_bins = hog_detector.nbins;
    float rad_per_bin = 3.14 / (float) n_bins;
    cv::Size win_size = hog_detector.winSize;
    cv::Size cell_size = hog_detector.cellSize;
    cv::Size block_size = hog_detector.blockSize;
    cv::Size block_stride = hog_detector.blockStride;

    // prepare data structure: 9 orientation / gradient strenghts for each cell
    int cells_in_x_dir = win_size.width / cell_size.width;
    int cells_in_y_dir = win_size.height / cell_size.height;
    int n_cells = cells_in_x_dir * cells_in_y_dir;
    int cells_per_block = (block_size.width / cell_size.width) * (block_size.height / cell_size.height);

    int blocks_in_x_dir = (win_size.width - block_size.width) / block_stride.width + 1;
    int blocks_in_y_dir = (win_size.height - block_size.height) / block_stride.height + 1;
    int n_blocks = blocks_in_x_dir * blocks_in_y_dir;

    float ***gradientStrengths = new float **[cells_in_y_dir];
    int **cellUpdateCounter = new int *[cells_in_y_dir];
    for (int y = 0; y < cells_in_y_dir; y++) {
        gradientStrengths[y] = new float *[cells_in_x_dir];
        cellUpdateCounter[y] = new int[cells_in_x_dir];
        for (int x = 0; x < cells_in_x_dir; x++) {
            gradientStrengths[y][x] = new float[n_bins];
            cellUpdateCounter[y][x] = 0;

            for (int bin = 0; bin < n_bins; bin++)
                gradientStrengths[y][x][bin] = 0.0;
        }
    }


//    compute gradient strengths per cell;
    int descriptorDataIdx = 0;


    for (int block_x = 0; block_x < blocks_in_x_dir; block_x++) {
        for (int block_y = 0; block_y < blocks_in_y_dir; block_y++) {
            int cell_start_x = block_x * block_stride.width / cell_size.width;
            int cell_start_y = block_y * block_stride.height / cell_size.height;

            for (int cell_id_x = cell_start_x;
                 cell_id_x < cell_start_x + block_size.width / cell_size.width; cell_id_x++)
                for (int cell_id_y = cell_start_y;
                     cell_id_y < cell_start_y + block_size.height / cell_size.height; cell_id_y++) {

                    for (int bin = 0; bin < n_bins; bin++) {
                        float val = feats.at(descriptorDataIdx++);
                        gradientStrengths[cell_id_y][cell_id_x][bin] += val;
                    }
                    cellUpdateCounter[cell_id_y][cell_id_x]++;
                }
        }
    }


    // compute average gradient strengths
    for (int celly = 0; celly < cells_in_y_dir; celly++) {
        for (int cellx = 0; cellx < cells_in_x_dir; cellx++) {

            float NrUpdatesForThisCell = (float) cellUpdateCounter[celly][cellx];

            // compute average gradient strenghts for each gradient bin direction
            for (int bin = 0; bin < n_bins; bin++) {
                gradientStrengths[celly][cellx][bin] /= NrUpdatesForThisCell;
            }
        }
    }


    for (int celly = 0; celly < cells_in_y_dir; celly++) {
        for (int cellx = 0; cellx < cells_in_x_dir; cellx++) {
            int drawX = cellx * cell_size.width;
            int drawY = celly * cell_size.height;

            int mx = drawX + cell_size.width / 2;
            int my = drawY + cell_size.height / 2;

            rectangle(visual_image,
                      cv::Point(drawX * scale_factor, drawY * scale_factor),
                      cv::Point((drawX + cell_size.width) * scale_factor,
                                (drawY + cell_size.height) * scale_factor),
                      CV_RGB(100, 100, 100),
                      1);

            for (int bin = 0; bin < n_bins; bin++) {
                float currentGradStrength = gradientStrengths[celly][cellx][bin];

                if (currentGradStrength == 0)
                    continue;

                float currRad = bin * rad_per_bin + rad_per_bin / 2;

                float dirVecX = cos(currRad);
                float dirVecY = sin(currRad);
                float maxVecLen = cell_size.width / 2;
                float scale = scale_factor / 5.0; // just a visual_imagealization scale,

                // compute line coordinates
                float x1 = mx - dirVecX * currentGradStrength * maxVecLen * scale;
                float y1 = my - dirVecY * currentGradStrength * maxVecLen * scale;
                float x2 = mx + dirVecX * currentGradStrength * maxVecLen * scale;
                float y2 = my + dirVecY * currentGradStrength * maxVecLen * scale;

                // draw gradient visual_imagealization
                line(visual_image,
                     cv::Point(x1 * scale_factor, y1 * scale_factor),
                     cv::Point(x2 * scale_factor, y2 * scale_factor),
                     CV_RGB(0, 0, 255),
                     1);

            }

        }
    }


    for (int y = 0; y < cells_in_y_dir; y++) {
        for (int x = 0; x < cells_in_x_dir; x++) {
            delete[] gradientStrengths[y][x];
        }
        delete[] gradientStrengths[y];
        delete[] cellUpdateCounter[y];
    }
    delete[] gradientStrengths;
    delete[] cellUpdateCounter;
    cv::imshow("HOG vis", visual_image);
    cv::waitKey(-1);
    cv::imwrite("hog_vis.jpg", visual_image);

}

