package ye.tian;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by ye.sensewhere on 2016/10/25.
 */
public class Demo2510 {

    private static final double betaDis = 1;
    private static final double betaSma = 0.2;
    private static final double betaOri = 2;
    private static final double eps = 1e-15;

    public static void main(String[] args) {

        long timeStamp = System.currentTimeMillis();

        FileIO fileIO = new FileIO();
        CoordinateTrans coordinateTrans = new CoordinateTrans();
        ShadowMatching shadowMatching = new ShadowMatching();

        String llaCsvPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/dataset/kdm6.lla.csv";
        ArrayList<double[]> llaArrayList = fileIO.llaLoader(llaCsvPath, ",");

        String satPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/dataset/kdm6.sat";
        ArrayList<ArrayList<double[]>> satArrList = fileIO.satelliteLoader(satPath, ",");

        String sensorsPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/dataset/kdm6.sen";
        ArrayList<double[]> sensorsArrList = fileIO.sensorsLoader(sensorsPath, ",");

        String bld2Dpath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/map/building_wangjing.MIF.out";
        String bldHgtPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/map/building_wangjing.MIF.height";
        ArrayList<BldModel> bldModels = fileIO.buildingModelLoader(bld2Dpath, bldHgtPath, ",");

//        System.out.println("Loading time: " + Long.toString(System.currentTimeMillis()-timeStamp) + " ms");

//        double[] llo = new double[]{39.9931,116.4684,0};
        double[] llo = llaArrayList.get(0);
//        ArrayList<ArrayList<double[][]>> bldLayoutArr = shadowMatching.getBldLayout(bldModels, llo, 200);

        ArrayList<double[]> states = new ArrayList<>();
        ArrayList<double[]> statesPre;
        ArrayList<double[]> llaOutArr = new ArrayList<>();
        ArrayList<ArrayList<double[][]>> bldLayoutArr = new ArrayList<>();
        double[] locOut = new double[]{0,0,0};
        double[] locInitial = new double[]{0,0,0};
        long gridShiftX, gridShiftY;

//        for (int t = 0, len = llaArrayList.size(); t < len; t++) {
        for (int t = 0, len = 8; t < len; t++) {
            double[] lla = llaArrayList.get(t);
            double sigma = llaArrayList.get(t)[2] * 10;
            ArrayList<double[]> satArr = satArrList.get(t);

//            double[] llo = lla;

            double[] locGps = coordinateTrans.llaToFlat(lla,llo,0,Math.PI/2);

//            double gridRadius = llaArrayList.get(t)[2] * 5;
            double gridRadius = 50;
            int gridStep = 2;

            if (t > 0) {

                double distPdr = 0.1 * (sensorsArrList.get(t)[2] - sensorsArrList.get(t - 1)[2]);
                double[] headingPdr = new double[]{Math.sin(sensorsArrList.get(t)[1]), -Math.cos(sensorsArrList.get(t)[1]), 0};

                statesPre = new ArrayList<>();
                gridShiftX = Math.round(locOut[0]/2) * 2;
                gridShiftY = Math.round(locOut[1]/2) * 2;
                for (double m = -gridRadius; m <= gridRadius; m += gridStep) {
                    for (double n = -gridRadius; n <= gridRadius; n += gridStep) {
                        double[] statePre = new double[]{locInitial[0] + gridShiftX + m, locInitial[1] + gridShiftY + n, 1, 1, 0, 0};
                        statesPre.add(statePre);
                    }
                }

//                bldLayoutArr = shadowMatching.getBldLayout(bldModels, llo, coordinateTrans.flatToLla(locOut,llo,0,Math.PI/2), 200);

                double sumStatesPre = 0;
                for (double[] statePre : statesPre) {
                    statePre[2] = 0;
                    for (double[] state : states) {
                        double distRoute = Math.sqrt((statePre[0] - state[0]) * (statePre[0] - state[0]) + (statePre[1] - state[1]) * (statePre[1] - state[1]));
                        double simDist = Math.exp(-Math.abs(distPdr - distRoute) / betaDis) / betaDis;
                        double simOri;
                        if (Arrays.equals(statePre,state)) {
                            System.out.println("Yes");
                            simOri = 1;
                        } else {
                            simOri = coordinateTrans.arrDotProd(headingPdr, coordinateTrans.arrSub(statePre, state)) / (coordinateTrans.arrNorm(headingPdr) * coordinateTrans.arrNorm(coordinateTrans.arrSub(statePre, state)));
                        }
                        simOri = 1 - 2 * Math.acos(simOri) / Math.PI;
                        if (simOri < 0.5) simOri = 0;
                        statePre[2] += state[5] * simDist * Math.pow(simOri, betaOri);
//                        if (statePre[2] < eps) statePre[2] = 0;
                    }
                    sumStatesPre += statePre[2];
                }

                for (double[] statePre : statesPre) {
                    statePre[2] /= sumStatesPre;
                }

                states = statesPre;

//                System.out.println("Transition time: " + Long.toString(System.currentTimeMillis()-timeStamp) + " ms");

            } else {

//                bldLayoutArr = shadowMatching.getBldLayout(bldModels, llo, llo, 200);

                for (double m = -gridRadius; m <= gridRadius; m += gridStep) {
                    for (double n = -gridRadius; n <= gridRadius; n += gridStep) {
                        double[] state = new double[]{locGps[0] + m, locGps[1] + n, 1, 1, 0, 0};
                        states.add(state);
                    }
                }

                locInitial = locGps;
//                System.out.println("Initialization time: " + Long.toString(System.currentTimeMillis()-timeStamp) + " ms");

//                for (double[] state : states) {
//                    state[2] = 1;
//                }
            }

            double[] matchScore = new double[states.size()];
            double matchScoreHighest = 0;
            double sumStatesGps = 0;
            for (int i = 0; i < states.size(); i++) {
//
                double[] state = states.get(i);
                state[3] = Math.exp(-0.5 * ((state[0] - locGps[0]) * (state[0] - locGps[0]) + (state[1] - locGps[1]) * (state[1] - locGps[1])) / sigma / sigma) / Math.sqrt(2 * Math.PI) / sigma;
                sumStatesGps += state[3];
//                if (state[3] < eps) state[3] = 0;
                bldLayoutArr =  shadowMatching.getBldLayout(bldModels, llo, coordinateTrans.flatToLla(new double[]{state[0],state[1],0},llo,0,Math.PI/2), 200);
                matchScore[i] = shadowMatching.getSmScore(new double[]{state[0], state[1], 0}, satArr, bldLayoutArr);
                if (matchScore[i] > matchScoreHighest) matchScoreHighest = matchScore[i];

//                System.out.println("Time " + Integer.toString(i) + " GPS and SM time: " + Long.toString(System.currentTimeMillis()-timeStamp) + " ms");
            }

            for (double[] state : states) {
                state[3] /= sumStatesGps;
            }

//            System.out.println("GPS and SM time: " + Long.toString(System.currentTimeMillis()-timeStamp) + " ms");

            double sumStates = 0;
            for (int i = 0; i < states.size(); i++) {
                states.get(i)[4] = Math.exp((matchScore[i] - matchScoreHighest) / betaSma) / betaSma;
//                if (states.get(i)[4] < eps) states.get(i)[4] = 0;
                states.get(i)[5] = states.get(i)[2] * states.get(i)[3] * states.get(i)[4];
//                if (states.get(i)[5] < eps) states.get(i)[5] = 0;
                sumStates += states.get(i)[5];
            }

            locOut = new double[]{0,0,0};
            for (double[] state : states) {
                state[5] /= sumStates;
                locOut[0] += state[0] * state[5];
                locOut[1] += state[1] * state[5];
            }

            llaOutArr.add(coordinateTrans.flatToLla(locOut,llo,0,Math.PI/2));

//            System.out.println("Output time: " + Long.toString(System.currentTimeMillis()-timeStamp) + " ms");

            BufferedWriter br = null;
            try {
                br = new BufferedWriter(new FileWriter("C:/Users/ye.sensewhere/Documents/new project on shadowing/D10192016/data/java/states_t23.csv", true));
                StringBuilder sb = new StringBuilder();
                for (double[] state : states) {
                    for (double item : state) {
                        sb.append(Double.toString(item));
                        sb.append(",");
                    }
                    sb.append("\n");
                }
                sb.append("\r\n");
                br.write(sb.toString());
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }

            System.out.println("Solved t = " + Integer.toString(t) + " ...");

        }

        BufferedWriter br = null;
        try {
            br = new BufferedWriter(new FileWriter("C:/Users/ye.sensewhere/Documents/new project on shadowing/D10192016/data/java/lla_out_t23.csv"));
            StringBuilder sb = new StringBuilder();
            for (double[] lla : llaOutArr) {
                for (double item : lla) {
                    sb.append(Double.toString(item));
                    sb.append(",");
                }
                sb.append("\n");
            }
            br.write(sb.toString());
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}