package ye.tian;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by ye.sensewhere on 2016/10/27.
 */
public class Demo2710 {

    private static final double betaDis = 0.1;
    private static final double betaSma = 0.2;
    private static final double betaOri = 2;
    private static final double eps = 1e-15;

    public static void main(String[] args) {

        long timeStamp = System.currentTimeMillis();

        FileIO fileIO = new FileIO();
        CoordinateTrans coordinateTrans = new CoordinateTrans();
        ShadowMatching shadowMatching = new ShadowMatching();

        String llaCsvPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/dataset/kdm4.lla.csv";
        ArrayList<double[]> llaArrayList = fileIO.llaLoader(llaCsvPath, ",");

        String satPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/dataset/kdm4.sat";
        ArrayList<ArrayList<double[]>> satArrList = fileIO.satelliteLoader(satPath, ",");

        String sensorsPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/dataset/kdm4.sen";
        ArrayList<double[]> sensorsArrList = fileIO.sensorsLoader(sensorsPath, ",");

        String bld2Dpath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/map/building_wangjing.MIF.out";
        String bldHgtPath = "C:/Users/ye.sensewhere/Documents/new project on shadowing/week39/map/building_wangjing.MIF.height";
        ArrayList<BldModel> bldModels = fileIO.buildingModelLoader(bld2Dpath, bldHgtPath, ",");

        System.out.println("Loading time: " + Long.toString(System.currentTimeMillis() - timeStamp) + " ms");

        boolean outDoor = false;
//        double[] llo = new double[]{39.9931,116.4684,0};
        double[] llo = llaArrayList.get(0);
        double[] locOut = new double[]{0, 0, 0};
        double[] locInitial = new double[]{0, 0, 0};
        long gridShiftX, gridShiftY;
//        ArrayList<ArrayList<double[][]>> bldLayoutArr = shadowMatching.getBldLayout(bldModels, llo, 200);

        ArrayList<double[]> states = new ArrayList<>();
        ArrayList<double[]> statesPre;
        ArrayList<double[]> llaOutArr = new ArrayList<>();

        for (int t = 0, len = llaArrayList.size(); t < len; t++) {

            double[] lla = llaArrayList.get(t);
            double sigma = llaArrayList.get(t)[2] * 10;
            ArrayList<double[]> satArr = satArrList.get(t);

            // TODO: 2016/10/27
            // GPS indoor/outdoor classification
            // GPS accuracy verification, may need 3D map layout information
            // By considering 3D layout information, to check if GPS accuracy indicator is accurate

            // AAR, current outdoor/indoor is just a temporary simple check by satellite info
            outDoor = (satArr.size() != 1 || satArr.get(0)[0] != 0);

//            double[] llo = lla;

            double[] locGps = coordinateTrans.llaToFlat(lla, llo, 0, Math.PI / 2);

//            double gridRadius = llaArrayList.get(t)[2] * 10;
            double gridRadius = 30;
            int gridStep = 2;

//            ArrayList<ArrayList<double[][]>> bldLayoutArr = shadowMatching.getBldLayout(bldModels, llo, llo, 200);
            ArrayList<ArrayList<double[][]>> bldLayoutArr = new ArrayList<>();

            if (t > 0) {

                double distPdr = 0.1 * (sensorsArrList.get(t)[2] - sensorsArrList.get(t-1)[2]);
                double[] headingPdr = new double[]{Math.sin(sensorsArrList.get(t)[1]), -Math.cos(sensorsArrList.get(t)[1]), 0};

                statesPre = new ArrayList<>();
                gridShiftX = Math.round(locOut[0] / gridStep) * gridStep;
                gridShiftY = Math.round(locOut[1] / gridStep) * gridStep;
                for (double m = -gridRadius; m <= gridRadius; m += gridStep) {
                    for (double n = -gridRadius; n <= gridRadius; n += gridStep) {
                        double[] statePre = new double[]{locInitial[0] + gridShiftX + m, locInitial[1] + gridShiftY + n, 1, 1, 1, 1};
                        if (MyMath.arrNorm(MyMath.arrSub(new double[]{statePre[0], statePre[1], 0}, new double[]{locInitial[0] + gridShiftX, locInitial[1] + gridShiftY, 0})) <= gridRadius)
                            statesPre.add(statePre);
                    }
                }

                bldLayoutArr = shadowMatching.getBldLayout(bldModels, llo, coordinateTrans.flatToLla(locOut, llo, 0, Math.PI / 2), 200);

                double sumStatesPre = 0;
                for (double[] statePre : statesPre) {
                    statePre[2] = 0;
                    for (double[] state : states) {
                        double distRoute = Math.sqrt((statePre[0] - state[0]) * (statePre[0] - state[0]) + (statePre[1] - state[1]) * (statePre[1] - state[1]));
                        double simDist = Math.exp(-Math.abs(distPdr - distRoute) / betaDis) / betaDis;
                        double simOri; // TODO: 2016/10/24 Check if Similarity is properly calculated
                        if (distRoute < eps) {
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

                System.out.println("Transition time: " + Long.toString(System.currentTimeMillis() - timeStamp) + " ms");

            } else {

                for (double m = -gridRadius; m <= gridRadius; m += gridStep) {
                    for (double n = -gridRadius; n <= gridRadius; n += gridStep) {
                        double[] state = new double[]{locGps[0] + m, locGps[1] + n, 1, 1, 1, 1};
                        if (MyMath.arrNorm(MyMath.arrSub(new double[]{state[0], state[1], 0}, locGps)) <= gridRadius)
                            states.add(state);
                    }
                }

                bldLayoutArr = shadowMatching.getBldLayout(bldModels, llo, llo, 200);

                locInitial = locGps;

                System.out.println("Initialization time: " + Long.toString(System.currentTimeMillis() - timeStamp) + " ms");
            }

            double[] matchScore = new double[states.size()];
            double matchScoreHighest = 0;
            double sumStatesGps = 0;
            for (int i = 0; i < states.size(); i++) {
//
                double[] state = states.get(i);
                if (outDoor) state[3] = Math.exp(-0.5 * ((state[0] - locGps[0]) * (state[0] - locGps[0]) + (state[1] - locGps[1]) * (state[1] - locGps[1])) / sigma / sigma) / Math.sqrt(2 * Math.PI) / sigma;
                sumStatesGps += state[3];
//                if (state[3] < eps) state[3] = 0;
//                bldLayoutArr =  shadowMatching.getBldLayout(bldModels, llo, coordinateTrans.flatToLla(new double[]{state[0],state[1],0},llo,0,Math.PI/2), 200);
                matchScore[i] = shadowMatching.getSmScore(new double[]{state[0], state[1], 0}, satArr, bldLayoutArr);
                if (matchScore[i] > matchScoreHighest) matchScoreHighest = matchScore[i];
            }

            double sumStatesSma = 0;
            for (int i = 0; i < states.size(); i++) {
                states.get(i)[4] = Math.exp((matchScore[i] - matchScoreHighest) / betaSma) / betaSma;
                sumStatesSma += states.get(i)[4];
            }

            for (double[] state : states) {
                state[3] /= sumStatesGps;
                state[4] /= sumStatesSma;
            }

            System.out.println("GPS and SM time: " + Long.toString(System.currentTimeMillis() - timeStamp) + " ms");

            double sumStates = 0;
            for (double[] state : states) {
                state[5] = state[2] * state[3] * state[4];
                sumStates += state[5];
            }

            locOut = new double[]{0, 0, 0};
            for (double[] state : states) {
                state[5] /= sumStates;
                locOut[0] += state[0] * state[5];
                locOut[1] += state[1] * state[5];
            }

            llaOutArr.add(coordinateTrans.flatToLla(locOut, llo, 0, Math.PI / 2));

            System.out.println("Combination time: " + Long.toString(System.currentTimeMillis() - timeStamp) + " ms");

            BufferedWriter br;
            try {
                br = new BufferedWriter(new FileWriter("C:/Users/ye.sensewhere/Documents/new project on shadowing/D10192016/data/java/states_kdm4.csv", true));
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

        BufferedWriter br;
        try {
            br = new BufferedWriter(new FileWriter("C:/Users/ye.sensewhere/Documents/new project on shadowing/D10192016/data/java/lla_out_kdm4.csv", true));
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

        System.out.println("Write to file time: " + Long.toString(System.currentTimeMillis() - timeStamp) + " ms");

    }
}