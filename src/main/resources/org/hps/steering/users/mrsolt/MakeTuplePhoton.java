<lcsim xmlns:xs="http://www.w3.org/2001/XMLSchema-instance" 
       xs:noNamespaceSchemaLocation="http://www.lcsim.org/schemas/lcsim/1.0/lcsim.xsd">
    <execute>
        <!--<driver name="RawTrackerHitSensorSetup"/>-->
    
        <driver name="PhotonTuple"/>
    </execute>    
    <drivers>    
        <driver name="RawTrackerHitSensorSetup" type="org.lcsim.recon.tracking.digitization.sisim.config.RawTrackerHitSensorSetup">
            <readoutCollections>SVTRawTrackerHits</readoutCollections>
        </driver>
        <driver name="PhotonTuple" type="org.hps.users.mrsolt.PhotonPositionDriver">
            <tupleFile>${outputFile}_photon.txt</tupleFile>
        </driver>
    </drivers>
</lcsim>

