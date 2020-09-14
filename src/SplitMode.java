import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class SplitMode {
    public static final int NUM_SQT_PER_TEMP = 20;

    public static void ClearTempDir(Path tempDirPath)
    {
        File tempDir = tempDirPath.toFile();
        for(File f: tempDir.listFiles())
        {
            if(Files.isSymbolicLink(f.toPath()))
            {
                f.delete();
            }
        }
        tempDir.delete();
    }

    public static Path CreateTempDirAndSoftLinkSqt(String parentPath, int numTemp, List<String> sqtList) throws IOException {
        String tempDirName = "temp_"+numTemp;
        Path tempPath = Paths.get(parentPath, tempDirName );
        if(tempPath.toFile().exists())
        {
            ClearTempDir(tempPath);
        }
        Files.createDirectory(tempPath);
        for(String sqt: sqtList)
        {
            Path sqtpath = Paths.get(parentPath, sqt);
            Path sqtTempPath = Paths.get(tempPath.toString(), sqt);
            Files.createSymbolicLink(sqtTempPath,sqtpath );
        }
        Path dtaParam = Paths.get(parentPath, "DTASelect.params");
        Path tempDTAParam = Paths.get(tempDirName, "DTASelect.params");
        Files.createSymbolicLink(tempDTAParam,dtaParam );

        Path sequestParam = Paths.get(parentPath, "sequest.params");
        Path tempSequestParam = Paths.get(tempDirName, "sequest.params");
        Files.createSymbolicLink(tempSequestParam,sequestParam );

        return tempPath;
    }

    public static List<String> FindAllSqt(String dirPath)
    {
        File dirF = new File(dirPath);
        List<String> sqtList = new LinkedList<>() ;
        String [] sqtArr = dirF.list((dir, name) ->name.endsWith(".sqt") );
        for(String s: sqtArr)
        {
            sqtList.add(s);
        }
        return sqtList;
    }

    public static void main(String [] args) throws IOException {
        String workingDir = System.getProperty("user.dir");
        List<String> sqtList= FindAllSqt(workingDir);
        int numTempDir = sqtList.size()/NUM_SQT_PER_TEMP;
        numTempDir = numTempDir%NUM_SQT_PER_TEMP == 0 ? numTempDir : numTempDir + 1;

        Collections.shuffle(sqtList);
        for(int i= 0; i<numTempDir; i++)
        {
            int end = (i+1)*NUM_SQT_PER_TEMP;
            end = end>=sqtList.size()? sqtList.size(): end;
            List<String> sqtSubList = sqtList.subList(i*NUM_SQT_PER_TEMP, end);
            CreateTempDirAndSoftLinkSqt(workingDir, i, sqtSubList);
        }
    }

}
