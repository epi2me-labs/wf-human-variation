
process pingMessage {
    label params.process_label
    cpus 1
    input:
        val message
        path json
    script:
        hostname = InetAddress.getLocalHost().getHostName()
        opsys = System.properties['os.name'].toLowerCase()
        disable = params.disable_ping ? '--disable' : ''
        meta = json.name != 'OPTIONAL_FILE' ? "--meta $json": ''
    """
    ping.py \
        --hostname $hostname \
        --opsys "$opsys" \
        --session $workflow.sessionId \
        --message $message \
        $meta $disable
    """
}


// send a start message
workflow start_ping {
    main:
       pingMessage("Started", Channel.fromPath("$projectDir/data/OPTIONAL_FILE"))
}

// send an end message
workflow end_ping {
    take:
        json
    main:
        pingMessage("Finished", json)
}
